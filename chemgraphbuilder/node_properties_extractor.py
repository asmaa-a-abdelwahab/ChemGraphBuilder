import os
import requests
import pandas as pd
import concurrent.futures
import xmltodict
from io import StringIO
import logging
import numpy as np
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor, as_completed
import math

logging.basicConfig(level=logging.INFO)

class NodePropertiesExtractor:
    """
    Extracts data from PubChem to build knowledge graphs in Neo4j, focusing on nodes representing chemical
    entities and their relationships. This class serves as a bridge between the PubChem database and Neo4j,
    allowing users to query chemical data and construct a graph-based representation of chemical compounds,
    their assays, related genes, and proteins.

    The primary functionality revolves around fetching detailed information about specified enzymes from PubChem,
    including assay data, gene properties, protein properties, and compound properties. It processes this data into
    a structured format suitable for knowledge graph construction, specifically tailored for use with Neo4j databases.

    Attributes:
        enzyme_list (list of str): Enzymes to query in the PubChem database.
        _base_url (str): Base URL for the PubChem API requests.
        _sep (str): Delimiter for parsing CSV data from PubChem.
        _enzyme_count (int): Number of enzymes in the enzyme_list, calculated at initialization.

    Parameters:
        enzyme_list (list of str): List of enzyme names for which assay data will be fetched from PubChem.
        base_url (str, optional): Base URL for PubChem API requests. Defaults to the assay target genesymbol endpoint.
        sep (str, optional): Separator used for parsing CSV data returned by PubChem. Defaults to ','.

    Usage Example:
        >>> enzyme_list = ['CYP2D6', 'CYP3A4']
        >>> extractor = NodePropertiesExtractor(enzyme_list)
        >>> df = extractor.run()
        This example initiates the extractor with a list of enzymes, fetches their data from PubChem,
        processes it, and potentially prepares it for knowledge graph construction in Neo4j.

    Note:
        To fully utilize this class, ensure you have network access to the PubChem API for data retrieval and
        a Neo4j database instance for knowledge graph construction. The class methods facilitate data extraction
        and processing, but integrating the output into Neo4j requires additional steps outside the scope of this class.
    """

    _REQUEST_TIMEOUT = 30  # in seconds
    _CONCURRENT_REQUEST_LIMIT = 2
    _RETRY_ATTEMPTS = 3  # number of times to retry a failed request

    def __init__(self, enzyme_list, base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol", sep=","):
        """
        Initializes a NodePropertiesExtractor instance, setting up the base URL for API requests, the separator for CSV parsing,
        and the list of enzymes to query from the PubChem database.

        Parameters:
            enzyme_list (list of str): A list of enzyme names for which to fetch assay data.
            base_url (str, optional): The base URL for PubChem API requests.
                                    Default is set to the assay target genesymbol endpoint.
            sep (str, optional): The delimiter to use for parsing CSV files returned by PubChem. Defaults to ','.

        Attributes:
            _base_url (str): Stores the base URL for API requests.
            _sep (str): Stores the delimiter for parsing CSV data.
            enzyme_list (list of str): Stores the list of enzyme names provided during initialization.
            _enzyme_count (int): The number of enzymes in the enzyme_list.
        """
        self._base_url = base_url
        self._sep = sep
        self.enzyme_list = enzyme_list
        self._enzyme_count = len(enzyme_list)

    def _make_request(self, url):
        """
        Sends an HTTP GET request to a specified URL with built-in retry logic. If the request fails,
        it retries the request up to a predefined number of attempts with exponential backoff to handle potential
        temporary network or server issues.

        The method attempts to gracefully handle server-side errors (HTTP 4XX/5XX responses) by raising an exception
        if the response status code indicates an error. For client-side errors (e.g., connectivity issues), it logs a warning
        and retries the request.

        Parameters:
            url (str): The complete URL to which the HTTP GET request is sent.

        Returns:
            requests.Response: The response object from the server if the request is successfully completed.

        Raises:
            requests.RequestException: If the request fails to complete successfully after the maximum number of retry attempts.
        """
        for attempt in range(self._RETRY_ATTEMPTS):
            try:
                response = requests.get(url, timeout=self._REQUEST_TIMEOUT)
                response.raise_for_status()  # Checks for HTTP errors
                return response
            except requests.RequestException as e:
                logging.warning(f"Attempt {attempt + 1} of {self._RETRY_ATTEMPTS} failed for URL: {url}. Error: {e}")
                if attempt + 1 == self._RETRY_ATTEMPTS:
                    raise  # All attempts failed; re-raise the last exception
                time.sleep(2 ** attempt)  # Exponential backoff


    def get_enzyme_assays(self, enzyme):
        """
        Fetches assay data for a specified enzyme from the PubChem database and returns it as a pandas DataFrame.

        This method constructs a URL to query the PubChem database for concise assay data related to the given enzyme.
        It processes the CSV response into a DataFrame, which includes various assay data points provided by PubChem.

        Parameters:
            enzyme (str): The name of the enzyme for which assay data is requested. This name is used in the API query.

        Returns:
            pd.DataFrame: A DataFrame containing the assay data fetched from PubChem for the specified enzyme. The DataFrame
                        includes columns based on the CSV response from PubChem, such as assay ID, results, and conditions.
                        Returns None if no data is available or if an error occurs during data fetching or processing.

        Raises:
            requests.RequestException: If an error occurs during the HTTP request to the PubChem API.
            pd.errors.EmptyDataError: If the response from PubChem contains no data.

        Example:
            >>> extractor = NodePropertiesExtractor(['enzyme'])
            >>> extractor.create_data_directories()
            >>> enzyme_assays_df = extractor.get_enzyme_assays('enzyme')
            >>> print(enzyme_assays_df.head())
        """
        assays_url = f"{self._base_url}/{enzyme}/concise/CSV"
        logging.info(f"Fetching assays for enzyme: {enzyme}")

        response = self._make_request(assays_url)

        assays_csv_string = response.text
        assays_csv_string_io = StringIO(assays_csv_string)
        try:
            assays_df = pd.read_csv(assays_csv_string_io, sep=self._sep, low_memory=False)
            logging.info(f"Assays DataFrame for enzyme {enzyme} has shape: {assays_df.shape}")
            return assays_df
        except pd.errors.EmptyDataError:
            logging.warning(f"No data available for enzyme {enzyme}.")
            return None


    def _process_enzymes(self, enzyme_list):
        """
        Iterates over a list of enzyme names, fetching assay data for each enzyme and aggregating the results into a list of DataFrames.

        This method calls `get_enzyme_assays` for each enzyme in the provided list, collecting the assay data (if available)
        into a list of pandas DataFrames. This list can then be used for further processing or analysis.

        Parameters:
            enzyme_list (list of str): A list containing the names of enzymes for which to fetch assay data.

        Returns:
            list of pd.DataFrame: A list containing a pandas DataFrame for each enzyme for which assay data was successfully
                                fetched and processed. Each DataFrame includes the assay data from PubChem for that enzyme.
                                If no data is available for an enzyme, it is omitted from the list.
        """
        df_list = [self.get_enzyme_assays(enzyme) for enzyme in enzyme_list]
        return [df for df in df_list if df is not None]


    def _concatenate_data(self, df_list):
        """
        Concatenates a list of pandas DataFrames into a single DataFrame. This method is useful for aggregating
        data fetched from multiple sources or APIs into a unified structure. If the list is empty, it returns None
        to indicate that no data was aggregated.

        Parameters:
            df_list (List[pd.DataFrame]): A list of pandas DataFrames to concatenate. These DataFrames should
                                        have the same structure (columns) to ensure proper concatenation.

        Returns:
            pd.DataFrame or None: A single concatenated DataFrame comprising all rows from the input DataFrames,
                                indexed continuously. Returns None if the input list is empty, indicating there
                                is no data to concatenate.
        """
        if df_list:
            return pd.concat(df_list, ignore_index=True)
        return None


    def run(self):
        """
        Orchestrates the process of fetching, filtering, and aggregating assay data from PubChem for a predefined list of enzymes.

        This method iteratively queries PubChem for assay data corresponding to each enzyme specified in the `enzyme_list` attribute during class initialization. It performs the following steps for each enzyme:
        1. Constructs a query URL and fetches assay data from PubChem.
        2. Filters the fetched data based on predefined criteria (e.g., containing specific substrings in the assay name).
        3. Aggregates the filtered data into a single pandas DataFrame.
        4. Identifies enzymes for which data could not be fetched or were excluded based on filtering criteria, logging their names.

        The final aggregated DataFrame, containing assay data for all successfully processed enzymes, is then saved to a CSV file. This method facilitates the extraction and preprocessing of chemical assay data for further analysis or integration into knowledge graphs.

        Note:
            - This method relies on the successful response from PubChem for each enzyme query.
            - Enzymes with no available data or failing to meet the filtering criteria are excluded from the final DataFrame.
            - The output CSV file is saved in the current working directory with the name 'Data/AllDataConnected.csv'.

        Returns:
            pd.DataFrame: A DataFrame containing the aggregated and filtered assay data for the specified enzymes.
                        Columns in the DataFrame correspond to the assay data fields returned by PubChem, subject to
                        the filtering criteria applied within this method.

        Raises:
            requests.RequestException: If there is an issue with fetching data from PubChem, such as a network problem or
                                    an invalid response.

        Example:
            Assuming `enzyme_list` was set to ['CYP2D6', 'CYP3A4'] during class initialization:

            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> result_df = extractor.run()
            >>> print(result_df.head())

            This will fetch and process assay data for 'CYP2D6' and 'CYP3A4', returning a DataFrame with the processed data.
        """

        # Initialize an empty list to store enzymes with successful responses
        enzymes_with_response = []

        # Keep a copy of the original list to identify removed enzymes later
        original_enzyme_list = self.enzyme_list.copy()

        for enzyme in self.enzyme_list:
            # Formulate the URL
            url = f"{self._base_url}/{enzyme}/concise/CSV"

            try:
                response = requests.get(url)
                # Check for a successful response (status code 200)
                if response.status_code == 200:
                    enzymes_with_response.append(enzyme)  # Keep the enzyme in the new list
            except requests.RequestException:
                # If there's an exception, skip adding the enzyme to the new list
                pass

        # Update the enzyme list with only the enzymes that had a successful response
        self.enzyme_list = enzymes_with_response

        # Identify and print the removed enzymes
        removed_enzymes = [enzyme for enzyme in original_enzyme_list if enzyme not in enzymes_with_response]
        if removed_enzymes:
            print("These enzymes were removed because their names aren't correct:", "[", ", ".join(removed_enzymes), "]")

        df_list = self._process_enzymes(self.enzyme_list)
        df = self._concatenate_data(df_list)
        substrings_to_filter = ['CYP', 'Cytochrome']
        pattern = '|'.join(substrings_to_filter)
        df = df[df['Assay Name'].str.contains(pattern, case=False, na=False)]
        df.to_csv('Data/AllDataConnected.csv', index=False)
        return df


    def _fetch_gene_details(self, gene_id):
        """
        Fetches gene details in parallel using the PubChem API.

        Args:
            gene_id (int): The gene ID for fetching details.

        Returns:
            tuple: Contains gene ID, symbol, taxonomy, taxonomy ID, and synonyms.
        """
        BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        url = f"{BASE_URL}/gene/geneid/{int(gene_id)}/summary/JSON"
        try:
            response = self._make_request(url)
            data = response.json()

            # Extracting the necessary details
            symbol = data['GeneSummaries']['GeneSummary'][0].get('Symbol', None)
            taxonomy = data['GeneSummaries']['GeneSummary'][0].get('Taxonomy', None)
            taxonomy_id = data['GeneSummaries']['GeneSummary'][0].get('TaxonomyID', None)
            synonyms = data['GeneSummaries']['GeneSummary'][0].get('Synonym', None)
            # print(type(synonyms))
            return gene_id, symbol, taxonomy, taxonomy_id, synonyms
        except Exception as e:
            logging.error(f"Error fetching details for gene_id {gene_id}: {e}")
            return gene_id, None, None, None, None


    def extract_gene_properties(self, main_data):
        """
        Extracts and processes gene properties from a given data source, specifically targeting genes
        relevant to the study (e.g., CYP enzymes) and records their details in a structured DataFrame.

        This method reads gene data from a CSV file specified by `main_data`, queries the PubChem database
        for additional properties of each unique gene ID found in the file, and compiles these properties
        into a new DataFrame. It focuses on fetching details like gene symbols, taxonomy, taxonomy IDs, and
        synonyms for each gene. The final DataFrame is filtered to include only genes of particular interest
        (e.g., certain CYP enzymes) and saved to a separate CSV file for further analysis or use.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which generated after running
                            `extractor.run()`.

        Returns:
            pd.DataFrame: A DataFrame containing the compiled gene properties, including GeneID, Symbol,
                        Taxonomy, Taxonomy ID, and Synonyms, filtered to include only specified genes
                        of interest. This DataFrame is also saved to 'Data/Nodes/Gene_Properties.csv'.

        Raises:
            Exception: If there's an issue reading the initial CSV file or fetching gene details from PubChem,
                    details of the exception are logged, and the method proceeds to process the next gene ID.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.run()
            >>> gene_properties_df = extractor.extract_gene_properties('Data/AllDataConnected.csv')
            >>> print(gene_properties_df.head())

            This would read gene IDs from 'Data/AllDataConnected.csv', fetch their properties from PubChem,
            and compile the details into a DataFrame, filtering for specified genes of interest and saving
            the results to 'Data/Nodes/Gene_Properties.csv'.

        Note:
            The method filters the resulting DataFrame to include only genes with symbols in the predefined
            enzyme_list. Adjust this list as necessary to match the focus of your study or application.
        """
        df = pd.read_csv(main_data)
        df_gene = pd.DataFrame(columns=['GeneID', 'Symbol', 'Taxonomy', 'Taxonomy ID', 'Synonyms'])

        unique_gene_ids = df['Target GeneID'].unique().tolist()

        gene_details = []

        for gene_id in unique_gene_ids:
            try:
                gene_id, symbol, taxonomy, taxonomy_id, synonyms = self._fetch_gene_details(gene_id)
                gene_details.append({
                    'GeneID': gene_id,
                    'Symbol': symbol,
                    'Taxonomy': taxonomy,
                    'Taxonomy ID': taxonomy_id,
                    'Synonyms': str(synonyms)
                })
            except Exception as exc:
                logging.error(f"Error occurred while processing gene_id {gene_id}: {exc}")
                gene_details.append({
                    'GeneID': gene_id,
                    'Symbol': None,
                    'Taxonomy': None,
                    'Taxonomy ID': None,
                    'Synonyms': None
                })

        # Now create the DataFrame from the list of dictionaries
        df_gene = pd.DataFrame(gene_details)
        n = self._enzyme_count
        gene_ids = df['Target GeneID'].value_counts().head(n).index.tolist()
        df_gene = df_gene[df_gene['GeneID'].isin([int(item) for item in gene_ids])]
        df_gene.to_csv('Data/Nodes/Gene_Properties.csv', sep=',', index=False)
        return df_gene


    def _fetch_assay_details(self, aid):
        """
        Fetches assay details from the PubChem API for a given assay ID.

        Args:
            aid (int): The assay ID to fetch details for.

        Returns:
            dict: A dictionary containing assay details like AID, SourceName, SourceID, Name, and Description.
                  Returns None if an error occurs during fetching or parsing.
        """
        BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        url = f"{BASE_URL}/assay/aid/{aid}/summary/XML"  # Constructing the API URL
        response = self._make_request(url)  # Making the API request
        xml_data = response.text  # Getting the response text

        try:
            # Parsing the XML response
            data_dict = xmltodict.parse(xml_data)
            properties = ['AID', 'SourceName', 'SourceID', 'Name', 'Description',
                          'Protocol', 'Comment', 'Method', 'Target', 'CIDCountAll',
                          'CIDCountActive', 'CIDCountInactive', 'CIDCountInconclusive',
                          'CIDCountUnspecified', 'CIDCountProbe']

            assay_data = {}
            # Extracting required properties from the parsed XML
            for prop in properties:
                assay_data[prop] = data_dict.get('AssaySummaries', {}).get('AssaySummary', {}).get(prop, None)
            return assay_data
        except Exception as e:
            logging.error(f"Error parsing XML for AID {aid}: {e}")
            return None


    def extract_assay_properties(self, main_data):
        """
        Extracts detailed properties of assays from PubChem for each unique assay ID found in the input data file.

        This method processes an input CSV file containing assay IDs (AID) and performs concurrent HTTP requests to
        fetch detailed assay properties from the PubChem database. The retrieved details include assay type, activity name,
        source name, source ID, name, and description. These properties are compiled into a new DataFrame, which is then
        saved to a CSV file for further analysis or use.

        The method employs a ThreadPoolExecutor to manage concurrent requests efficiently, improving the performance
        when dealing with a large number of assay IDs. Errors encountered during data fetching are logged, and the
        process continues with the next assay ID, ensuring the method's robustness.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which generated after running
                            `extractor.run()`.

        Returns:
            pd.DataFrame: A DataFrame containing the fetched assay properties, including columns for AID, Assay Type,
                        Activity Name, SourceName, SourceID, Name, and Description. This DataFrame is saved to
                        'Data/Nodes/Assay_Properties.csv' in the current working directory.

        Raises:
            ValueError: If the input CSV file is empty or does not contain the 'AID' column.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.run()
            >>> assay_properties_df = extractor.extract_assay_properties('Data/AllDataConnected.csv')
            >>> print(assay_properties_df.head())

            This example reads assay IDs from 'Data/AllDataConnected.csv', queries PubChem for their detailed properties,
            and compiles the results into a DataFrame, which is also saved to 'Data/Nodes/Assay_Properties.csv'.

        Note:
            This method requires network access to the PubChem API and assumes the availability of a valid 'AID' column
            in the input CSV file. Ensure the input file path is correct and accessible to avoid errors during processing.
        """

        df = pd.read_csv(main_data)

        # Check if the DataFrame is valid
        if df.empty or 'AID' not in df.columns:
            logging.error("DataFrame is empty or does not contain 'AID' column.")
            return pd.DataFrame()

        unique_aids = df['AID'].unique().tolist()  # Extracting unique assay IDs

        columns = ['AID', 'Assay Type', 'Activity Name', 'SourceName',
                   'SourceID', 'Name', 'Description', 'Protocol',
                   'Comment', 'Method', 'Target', 'CIDCountAll',
                   'CIDCountActive', 'CIDCountInactive', 'CIDCountInconclusive', 
                   'CIDCountUnspecified', 'CIDCountProbe']
        assay_df = pd.DataFrame(columns=columns)  # Initializing a DataFrame to store assay properties

        # Using ThreadPoolExecutor for concurrent fetching of assay details
        with concurrent.futures.ThreadPoolExecutor(max_workers=self._CONCURRENT_REQUEST_LIMIT) as executor:
            future_to_aid = {executor.submit(self._fetch_assay_details, aid): aid for aid in unique_aids}

            # Iterating over completed futures
            for future in concurrent.futures.as_completed(future_to_aid):
                aid = future_to_aid[future]
                try:
                    assay_data = future.result()  # Fetching the result from the future
                    if assay_data:
                        # Preparing a new row with the fetched data
                        new_row = {
                            'AID': aid,
                            'Assay Type': df.loc[df['AID'] == aid, 'Assay Type'].iloc[0],
                            'Activity Name': df.loc[df['AID'] == aid, 'Activity Name'].iloc[0],
                            **assay_data
                        }
                        # Adding the new row to the DataFrame
                        assay_df = pd.concat([assay_df, pd.DataFrame([new_row])], ignore_index=True)
                except Exception as exc:
                    # Logging any errors encountered during the fetch
                    logging.error(f"Error occurred while processing AID {aid}: {exc}")

        # Saving the updated DataFrame to a CSV file
        assay_df.to_csv('Data/Nodes/Assay_Properties.csv', sep=',', index=False)
        return assay_df


    def extract_protein_properties(self, main_data):
        """
        Extracts and compiles protein properties from the NCBI protein database based on accession numbers.

        Given a CSV file specified by `main_data`, this method reads protein accession numbers and performs web
        scraping on the NCBI protein database pages to extract protein titles. The method constructs a URL for
        each accession number, sends a request to retrieve the page content, and parses the HTML to find the
        protein title. The extracted titles, along with their corresponding accession numbers and URLs, are
        compiled into a DataFrame. This DataFrame is saved to a CSV file, providing a structured summary of
        protein properties for further analysis or use.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which generated after running
                            `extractor.run()`.

        Returns:
            pd.DataFrame: A DataFrame with columns 'RefSeq Accession', 'URL', and 'Description', where
                        'Description' contains the title of the protein extracted from its NCBI page.
                        This DataFrame is saved to 'Data/Nodes/Protein_Properties.csv' in the current
                        working directory.

        Raises:
            Exception: If there's an issue reading the initial CSV file or querying the NCBI database,
                    details of the exception are logged. The method continues processing the next
                    accession number, ensuring robustness against individual failures.

        Example:
            Assuming 'protein_data.csv' contains a column 'Target Accession' with accession numbers:

            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.run() # you need to run this only once
            >>> protein_properties_df = extractor.extract_protein_properties('Data/AllDataConnected.csv')
            >>> print(protein_properties_df.head())

            This would read accession numbers from 'Data/AllDataConnected.csv', scrape their titles from the
            NCBI protein database, and compile the results into a DataFrame, which is also saved to
            'Data/Nodes/Protein_Properties.csv'.

        Note:
            This method requires internet access to query the NCBI protein database. Ensure the input file
            path is correct and accessible to avoid errors during processing. Web scraping is dependent on
            the structure of the web page; changes to the NCBI protein database pages may require updates
            to the scraping logic.
        """

        # Initialize a list to store the extracted data
        data = []

        n = self._enzyme_count
        df = pd.read_csv(main_data)
        gene_ids = df['Target GeneID'].value_counts().head(n).index.tolist()
        df = df[df['Target GeneID'].isin([int(item) for item in gene_ids])]
        Accessions = df['Target Accession'].unique().tolist()
        # Iterate over each protein accession number in the DataFrame
        for accession in Accessions:
            # Construct the URL to query the NCBI protein database
            url = f"https://www.ncbi.nlm.nih.gov/protein/{accession}"

            try:
                # Send an HTTP request to the URL
                response = requests.get(url)

                # Parse the HTML content of the response
                soup = BeautifulSoup(response.text, 'html.parser')

                # Extract the title from the parsed HTML
                title = soup.title.string if soup.title else 'Title Not Found'

                # Append the extracted data to the list
                data.append({'RefSeq Accession': accession, 'URL': url, 'Description': title})
            except Exception as e:
                # In case of an error, log the error message
                data.append({'RefSeq Accession': accession, 'URL': url, 'Description': f'Error: {e}'})

        # Convert the list of data into a DataFrame
        protein_df = pd.DataFrame(data)

        # Save the DataFrame to a CSV file
        protein_df.to_csv('Data/Nodes/Protein_Properties.csv', sep=',', index=False)

        # Return the DataFrame
        return protein_df


    def fetch_data(self, cid):
        """
        Retrieves detailed chemical compound properties for a specified Compound ID (CID) from the PubChem database.

        This method constructs a query URL to fetch a wide range of properties for the given CID from PubChem, including
        molecular formula, molecular weight, canonical and isomeric SMILES, InChI codes, physicochemical properties, and
        more. If the CID is valid and data is available, it returns a pandas DataFrame containing these properties. This
        method also generates a URL to retrieve the structure image of the compound as a 2D PNG image, adding it as a
        column in the DataFrame. In cases where the CID is NaN or an error occurs during data retrieval, an empty DataFrame
        is returned.

        Parameters:
            cid (int or float): The Compound ID for which to fetch data. Can be an integer or NaN.

        Returns:
            pd.DataFrame: A DataFrame containing the fetched properties for the given CID. The DataFrame includes
                        columns for each property fetched from PubChem, along with a 'StructureImage2DURL' column
                        containing the URL to the compound's structure image. Returns an empty DataFrame if the CID
                        is NaN or if any error occurs during the fetch operation.

        Raises:
            Exception: Logs an error message if the request to PubChem fails or if the response cannot be processed
                    into a DataFrame.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> compound_data_df = extractor.fetch_data(2244)
            >>> print(compound_data_df.head())

            This example fetches the properties for the compound with CID 2244 from PubChem and prints the first few rows
            of the resulting DataFrame.

        Note:
            This method requires an active internet connection to access the PubChem database. Ensure that the CID provided
            is valid and not NaN to avoid fetching errors. The structure and availability of data fields are subject to the
            current state of the PubChem database and may vary.
        """
        if pd.isna(cid):
            return pd.DataFrame()  # Return an empty DataFrame for NaN CIDs

        cid = int(cid)  # Convert CID to integer
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
               "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,"
               "InChIKey,IUPACName,Title,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,"
               "Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,"
               "IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,"
               "BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,"
               "PatentCount,PatentFamilyCount,LiteratureCount,Volume3D,XStericQuadrupole3D,"
               "YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,"
               "FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,"
               "FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D,"
               "Fingerprint2D/CSV")
        try:
            response = requests.get(url)
            response.raise_for_status()
            compound_data = pd.read_csv(StringIO(response.text), sep=',', low_memory=False)
            compound_data['StructureImage2DURL'] = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
            return compound_data
        except Exception as e:
            logging.error(f"Error processing CID {cid}: {e}")
            return pd.DataFrame()  # Return an empty DataFrame in case of error


    def extract_compound_properties(self, main_data):
        """
        Extracts and aggregates compound properties from PubChem for a list of compounds associated with specific genes.

        This method processes a CSV file specified by `main_data`, which contains gene identifiers and their associated
        compound IDs (CIDs). It selects compounds related to the top `n` most frequently occurring genes in the dataset,
        where `n` is determined by the instance's `_enzyme_count` attribute. The method then fetches detailed compound
        properties from PubChem in chunks, using concurrent requests to improve efficiency and manage the load on the
        PubChem API. The fetched compound properties are aggregated into a single DataFrame and saved to multiple CSV files,
        one for each chunk of compound IDs processed.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which generated after running
                            `extractor.run()`.

        Side Effects:
            - Saves the aggregated compound properties to CSV files in the current working directory. The files are named
            'Data/Nodes/Compound_Properties/Chunk_{i}.csv', where `{i}` is the chunk index.

        Returns:
            None: This method does not return a value. Instead, it saves the fetched compound data directly to CSV files.

        Raises:
            Exception: Logs an error and continues processing the next CID if an error occurs while fetching data for a specific CID.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.extract_compound_properties('Data/AllDataConnected.csv')
            This will read 'Data/AllDataConnected.csv', filter for compounds associated with the top n genes, fetch their properties
            from PubChem, and save the results into multiple CSV files for each chunk of compounds processed.

        Note:
            - Ensure that the 'main_data' CSV file exists and is accessible at the specified path.
            - The method automatically handles NaN values in the 'CID' column and excludes them from processing.
            - The `enzyme_count` attribute determines the number of top genes for which compound properties will be fetched.
            - Internet access is required to fetch compound data from the PubChem API.
            - The method employs a `ThreadPoolExecutor` with a configurable number of workers (default is len(enzyme_list)) to parallelize
            requests, which can be adjusted based on system capabilities and API rate limits.
        """

        n = self._enzyme_count
        df = pd.read_csv(main_data)
        gene_ids = df['Target GeneID'].value_counts().head(n).index.tolist()
        df = df[df['Target GeneID'].isin([int(item) for item in gene_ids])]
        df = df.dropna(subset=['CID'])
        IDs = df['CID'].unique().tolist()

        # Define chunk size and calculate number of chunks
        chunk_size = 10000
        num_chunks = math.ceil(len(IDs) / chunk_size)

        for i in range(num_chunks):
            # Calculate start and end indices for each chunk
            start_index = i * chunk_size
            end_index = start_index + chunk_size

            # Extract chunk of CIDs
            chunk_cids = IDs[start_index:end_index]
            # chunk_cids = [x for x in chunk_cids if not np.isnan(x)]

            # Use ThreadPoolExecutor to parallelize requests for the chunk
            with ThreadPoolExecutor(max_workers=5) as executor:
                future_to_cid = {executor.submit(self.fetch_data, cid): cid for cid in chunk_cids}
                results = []

                for future in as_completed(future_to_cid):
                    cid = future_to_cid[future]
                    try:
                        data = future.result()
                        results.append(data)
                    except Exception as e:
                        logging.error(f"Error processing CID {cid}: {e}")

            # Concatenate results for the current chunk
            chunk_df = pd.concat(results, ignore_index=True)

            # Save the concatenated DataFrame to a CSV file for the chunk
            chunk_df.to_csv(f'Data/Nodes/Compound_Properties/Chunk_{i}.csv', sep=',', index=False)

        return
