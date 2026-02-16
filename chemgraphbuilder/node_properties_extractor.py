"""
This module defines the `NodePropertiesExtractor` class, responsible for
extracting data from the PubChem database to build knowledge graphs in Neo4j.
The class focuses on nodes representing chemical entities and their relationships,
allowing users to query chemical data and construct a graph-based representation
of chemical compounds, their assays, related genes, and proteins.

The primary functionality revolves around fetching detailed information about
specified enzymes from PubChem, including assay data, gene properties, protein
properties, and compound properties. It processes this data into a structured
format suitable for knowledge graph construction, specifically tailored for use
with Neo4j databases.

Classes:
    - NodePropertiesExtractor: A class to extract data from PubChem to build
      knowledge graphs in Neo4j.

Usage Example:
    >>> enzyme_list = ['CYP2D6', 'CYP3A4']
    >>> extractor = NodePropertiesExtractor(enzyme_list)
    >>> df = extractor.run()
    This example initiates the extractor with a list of enzymes, fetches their
    data from PubChem, processes it, and potentially prepares it for knowledge
    graph construction in Neo4j.

Note:
    To fully utilize this class, ensure you have network access to the PubChem
    API for data retrieval and a Neo4j database instance for knowledge graph
    construction. The class methods facilitate data extraction and processing,
    but integrating the output into Neo4j requires additional steps outside the
    scope of this class.
"""

import datetime
import random
import time
import math
import logging
from io import StringIO
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Union
from bs4 import BeautifulSoup
import requests
import pandas as pd
import xmltodict

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class NodePropertiesExtractor:
    """
    Extracts data from PubChem to build knowledge graphs in Neo4j,
    focusing on nodes representing chemical entities and their relationships.
    This class serves as a bridge between the PubChem database and Neo4j,
    allowing users to query chemical data and construct a graph-based
    representation of chemical compounds, their assays, related genes, and proteins.

    The primary functionality revolves around fetching detailed information
    about specified enzymes from PubChem, including assay data, gene properties,
    protein properties, and compound properties. It processes this data into
    a structured format suitable for knowledge graph construction, specifically
    tailored for use with Neo4j databases.

    Attributes:
        enzyme_list (list of str): Enzymes to query in the PubChem database.
        _base_url (str): Base URL for the PubChem API requests.
        _sep (str): Delimiter for parsing CSV data from PubChem.
        _enzyme_count (int): Number of enzymes in the enzyme_list, calculated at
        initialization.

    Parameters:
        enzyme_list (list of str): List of enzyme names for which assay data
        will be fetched from PubChem.
        base_url (str, optional): Base URL for PubChem API requests. Defaults to
        the assay target genesymbol endpoint.
        sep (str, optional): Separator used for parsing CSV data returned
        by PubChem. Defaults to ','.

    Usage Example:
        >>> enzyme_list = ['CYP2D6', 'CYP3A4']
        >>> extractor = NodePropertiesExtractor(enzyme_list)
        >>> df = extractor.run()
        This example initiates the extractor with a list of enzymes, fetches
        their data from PubChem, processes it, and potentially prepares it for
        knowledge graph construction in Neo4j.

    Note:
        To fully utilize this class, ensure you have network access to the
        PubChem API for data retrieval and a Neo4j database instance for
        knowledge graph construction. The class methods facilitate data extraction
        and processing, but integrating the output into Neo4j requires additional
        steps outside the scope of this class.
    """

    _REQUEST_TIMEOUT = 300  # in seconds
    _CONCURRENT_REQUEST_LIMIT = 2
    _RETRY_ATTEMPTS = 3  # number of times to retry a failed request

    def __init__(
        self,
        enzyme_list,
        base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol",
        sep=",",
    ):
        """
        Initializes a NodePropertiesExtractor instance, setting up the base URL
        for API requests, the separator for CSV parsing, and the list of enzymes
        to query from the PubChem database.

        Parameters:
            enzyme_list (list of str): A list of enzyme names for which to fetch
            assay data.
            base_url (str, optional): The base URL for PubChem API requests.
            Default is set to the assay target genesymbol endpoint.
            sep (str, optional): The delimiter to use for parsing CSV files
            returned by PubChem. Defaults to ','.

        Attributes:
            _base_url (str): Stores the base URL for API requests.
            _sep (str): Stores the delimiter for parsing CSV data.
            enzyme_list (list of str): Stores the list of enzyme names provided
            during initialization.
            _enzyme_count (int): The number of enzymes in the enzyme_list.
        """
        self._base_url = base_url
        self._sep = sep
        self.enzyme_list = enzyme_list
        self._enzyme_count = len(enzyme_list)

    def _make_request(self, url):
        for attempt in range(self._RETRY_ATTEMPTS):
            try:
                response = requests.get(url, timeout=self._REQUEST_TIMEOUT)
                response.raise_for_status()
                return response
            except requests.HTTPError as e:
                status = getattr(e.response, "status_code", None)
                # For PubChem 503/504, treat as transient
                if status in (503, 504):
                    logging.warning(
                        "Attempt %s of %s failed for URL: %s (HTTP %s). Error: %s",
                        attempt + 1, self._RETRY_ATTEMPTS, url, status, e
                    )
                else:
                    # Non-transient errors: re-raise immediately
                    raise

            except requests.RequestException as e:
                logging.warning(
                    "Attempt %s of %s failed for URL: %s. Error: %s",
                    attempt + 1, self._RETRY_ATTEMPTS, url, e
                )

            # If we got here, there was some error
            if attempt + 1 == self._RETRY_ATTEMPTS:
                raise

            # Jittered backoff to be nicer to the API
            # sleep_seconds = (2 ** attempt) + random.uniform(0, 0.5)
            time.sleep(0.2)

    def get_enzyme_assays(self, enzyme):
        """
        Fetches assay data for a specified enzyme from the PubChem database and
        returns it as a pandas DataFrame.

        This method constructs a URL to query the PubChem database for concise
        assay data related to the given enzyme. It processes the CSV response
        into a DataFrame, which includes various assay data points provided by PubChem.

        Parameters:
            enzyme (str): The name of the enzyme for which assay data is
            requested. This name is used in the API query.

        Returns:
            pd.DataFrame: A DataFrame containing the assay data fetched from
            PubChem for the specified enzyme. The DataFrame includes columns
            based on the CSV response from PubChem, such as assay ID, results,
            and conditions. Returns None if no data is available or if an error
            occurs during data fetching or processing.

        Raises:
            requests.RequestException: If an error occurs during the HTTP
            request to the PubChem API.
            pd.errors.EmptyDataError: If the response from PubChem contains no data.

        Example:
            >>> extractor = NodePropertiesExtractor(['enzyme'])
            >>> enzyme_assays_df = extractor.get_enzyme_assays('enzyme')
            >>> print(enzyme_assays_df.head())
        """
        assays_url = f"{self._base_url}/{enzyme.lower()}/concise/CSV"
        logging.info(f"Fetching assays for enzyme: {enzyme}")

        response = self._make_request(assays_url)

        assays_csv_string = response.text
        assays_csv_string_io = StringIO(assays_csv_string)
        try:
            assays_df = pd.read_csv(
                assays_csv_string_io, sep=self._sep, low_memory=False
            )
            logging.info(
                "Assays DataFrame for enzyme %s has shape: %s", enzyme, assays_df.shape
            )
            return assays_df
        except pd.errors.EmptyDataError:
            logging.warning("No data available for enzyme %s.", enzyme)
            return None

    def _process_enzymes(self, enzyme_list):
        """
        Iterates over a list of enzyme names, fetching assay data for each enzyme
        and aggregating the results into a list of DataFrames.

        This method calls `get_enzyme_assays` for each enzyme in the provided
        list, collecting the assay data (if available) into a list of pandas
        DataFrames. This list can then be used for further processing or analysis.

        Parameters:
            enzyme_list (list of str): A list containing the names of enzymes
            for which to fetch assay data.

        Returns:
            list of pd.DataFrame: A list containing a pandas DataFrame for each
            enzyme for which assay data was successfully fetched and processed.
            Each DataFrame includes the assay data from PubChem for that enzyme.
            If no data is available for an enzyme, it is omitted from the list.
        """
        df_list = [self.get_enzyme_assays(enzyme) for enzyme in enzyme_list]
        return [df for df in df_list if df is not None]

    def _concatenate_data(self, df_list):
        """
        Concatenates a list of pandas DataFrames into a single DataFrame.
        This method is useful for aggregating
        data fetched from multiple sources or APIs into a unified structure.
        If the list is empty, it returns None to indicate that no data was
        aggregated.

        Parameters:
            df_list (List[pd.DataFrame]): A list of pandas DataFrames to
            concatenate. These DataFrames should have the same structure
            (columns) to ensure proper concatenation.

        Returns:
            pd.DataFrame or None: A single concatenated DataFrame comprising all
            rows from the input DataFrames, indexed continuously. Returns None
            if the input list is empty, indicating there is no data to concatenate.
        """
        if df_list:
            return pd.concat(df_list, ignore_index=True)
        return None

    def run(self):
        """
        Orchestrates the process of fetching, filtering, and aggregating assay
        data from PubChem for a predefined list of enzymes.

        This method iteratively queries PubChem for assay data corresponding
        to each enzyme specified in the `enzyme_list` attribute during class
        initialization. It performs the following steps for each enzyme:
        1. Constructs a query URL and fetches assay data from PubChem.
        2. Filters the fetched data based on predefined criteria
        (e.g., containing specific substrings in the assay name).
        3. Aggregates the filtered data into a single pandas DataFrame.
        4. Identifies enzymes for which data could not be fetched or were
        excluded based on filtering criteria, logging their names.

        The final aggregated DataFrame, containing assay data for all successfully
        processed enzymes, is then saved to a CSV file. This method facilitates
        the extraction and preprocessing of chemical assay data for further
        analysis or integration into knowledge graphs.

        Note:
            - This method relies on the successful response from PubChem
            for each enzyme query.
            - Enzymes with no available data or failing to meet the filtering
            criteria are excluded from the final DataFrame.
            - The output CSV file is saved in the current working directory
            with the name 'Data/AllDataConnected.csv'.

        Returns:
            pd.DataFrame: A DataFrame containing the aggregated and filtered
            assay data for the specified enzymes. Columns in the DataFrame
            correspond to the assay data fields returned by PubChem, subject to
            the filtering criteria applied within this method.

        Raises:
            requests.RequestException: If there is an issue with fetching data
            from PubChem, such as a network problem or an invalid response.

        Example:
            Assuming `enzyme_list` was set to ['CYP2D6', 'CYP3A4'] during
            class initialization:

            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> result_df = extractor.run()
            >>> print(result_df.head())

            This will fetch and process assay data for 'CYP2D6' and 'CYP3A4',
            returning a DataFrame with the processed data.
        """

        # Initialize an empty list to store enzymes with successful responses
        enzymes_with_response = []

        # Keep a copy of the original list to identify removed enzymes later
        original_enzyme_list = self.enzyme_list.copy()

        for enzyme in self.enzyme_list:
            url = f"{self._base_url}/{enzyme}/concise/CSV"
            try:
                response = self._make_request(url)
                if response.status_code == 200:
                    enzymes_with_response.append(enzyme)
            except requests.RequestException:
                # If there's an exception, skip adding the enzyme to the new list
                pass

        # Update the enzyme list with only the enzymes that had a successful response
        self.enzyme_list = enzymes_with_response

        # Identify and print the removed enzymes
        removed_enzymes = [
            enzyme
            for enzyme in original_enzyme_list
            if enzyme not in enzymes_with_response
        ]
        if removed_enzymes:
            logging.info(
                "These enzymes were removed because their names aren't correct: %s",
                ", ".join(removed_enzymes),
            )

        df_list = self._process_enzymes(self.enzyme_list)
        df = self._concatenate_data(df_list)
        substrings_to_filter = ["CYP", "Cytochrome"]
        pattern = "|".join(substrings_to_filter)
        df = df[df["Assay Name"].str.contains(pattern, case=False, na=False)]
        retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"
        df["SourceDatabase"] = "PubChem"
        df["SourceEndpoint"] = self._base_url
        df["RetrievedAt"] = retrieval_time
        df.to_csv("Data/AllDataConnected.csv", index=False)
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
            symbol = data["GeneSummaries"]["GeneSummary"][0].get("Symbol", None)
            taxonomy = data["GeneSummaries"]["GeneSummary"][0].get("Taxonomy", None)
            taxonomy_id = data["GeneSummaries"]["GeneSummary"][0].get(
                "TaxonomyID", None
            )
            synonyms = data["GeneSummaries"]["GeneSummary"][0].get("Synonym", None)
            # print(type(synonyms))
            return gene_id, symbol, taxonomy, taxonomy_id, synonyms
        except Exception as e:
            logging.error(f"Error fetching details for gene_id {gene_id}: {e}")
            return gene_id, None, None, None, None

    def extract_gene_properties(self, main_data):
        """
        Extracts and processes gene properties from a given data source,
        specifically targeting genes relevant to the study (e.g., CYP enzymes)
        and records their details in a structured DataFrame.

        This method reads gene data from a CSV file specified by `main_data`,
        queries the PubChem database for additional properties of each unique
        gene ID found in the file, and compiles these properties into a new
        DataFrame. It focuses on fetching details like gene symbols, taxonomy,
        taxonomy IDs, and synonyms for each gene. The final DataFrame is filtered
        to include only genes of particular interest (e.g., certain CYP enzymes)
        and saved to a separate CSV file for further analysis or use.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which
            generated after running `extractor.run()`.

        Returns:
            pd.DataFrame: A DataFrame containing the compiled gene properties,
            including GeneID, Symbol, Taxonomy, Taxonomy ID, and Synonyms,
            filtered to include only specified genes of interest. This DataFrame
            is also saved to 'Data/Nodes/Gene_Properties.csv'.

        Raises:
            Exception: If there's an issue reading the initial CSV file or
            fetching gene details from PubChem, details of the exception are
            logged, and the method proceeds to process the next gene ID.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.run()
            >>> gene_properties_df = extractor.extract_gene_properties('Data/AllDataConnected.csv')
            >>> print(gene_properties_df.head())

            This would read gene IDs from 'Data/AllDataConnected.csv', fetch
            their properties from PubChem, and compile the details into a
            DataFrame, filtering for specified genes of interest and saving
            the results to 'Data/Nodes/Gene_Properties.csv'.

        Note:
            The method filters the resulting DataFrame to include only genes with
            symbols in the predefined enzyme_list. Adjust this list as necessary
            to match the focus of your study or application.
        """
        df = pd.read_csv(main_data)
        df_gene = pd.DataFrame(
            columns=["GeneID", "Symbol", "Taxonomy", "Taxonomy ID", "Synonyms"]
        )

        unique_gene_ids = df["Target GeneID"].unique().tolist()

        gene_details = []

        for gene_id in unique_gene_ids:
            try:
                gene_id, symbol, taxonomy, taxonomy_id, synonyms = (
                    self._fetch_gene_details(gene_id)
                )
                gene_details.append(
                    {
                        "GeneID": gene_id,
                        "Symbol": symbol,
                        "Taxonomy": taxonomy,
                        "Taxonomy ID": taxonomy_id,
                        "Synonyms": str(synonyms),
                    }
                )
            except Exception as exc:
                logging.error(
                    "Error occurred while processing gene_id %s: %s", gene_id, exc
                )
                gene_details.append(
                    {
                        "GeneID": gene_id,
                        "Symbol": None,
                        "Taxonomy": None,
                        "Taxonomy ID": None,
                        "Synonyms": None,
                    }
                )

        # Now create the DataFrame from the list of dictionaries
        df_gene = pd.DataFrame(gene_details)
        n = self._enzyme_count
        gene_ids = df["Target GeneID"].value_counts().head(n).index.tolist()
        df_gene = df_gene[df_gene["GeneID"].isin([int(item) for item in gene_ids])]
        # after df_gene is created and filtered
        retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"
        df_gene["SourceDatabase"] = "PubChem"
        df_gene["SourceEndpoint"] = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/summary/JSON"
        df_gene["RetrievedAt"] = retrieval_time

        df_gene.to_csv("Data/Nodes/Gene_Properties.csv", sep=",", index=False)

        return df_gene

    def _fetch_assay_details(self, aid):
        """
        Fetches assay details from the PubChem API for a given assay ID.

        Args:
            aid (int): The assay ID to fetch details for.

        Returns:
            dict: A dictionary containing assay details like AID, SourceName,
            SourceID, Name, and Description. Returns None if an error occurs
            during fetching or parsing.
        """
        BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        url = f"{BASE_URL}/assay/aid/{aid}/summary/XML"  # Constructing the API URL
        response = self._make_request(url)  # Making the API request
        xml_data = response.text  # Getting the response text

        try:
            # Parsing the XML response
            data_dict = xmltodict.parse(xml_data)
            properties = [
                "AID",
                "SourceName",
                "SourceID",
                "Name",
                "Description",
                "Protocol",
                "Comment",
                "Method",
                "Target",
                "CIDCountAll",
                "CIDCountActive",
                "CIDCountInactive",
                "CIDCountInconclusive",
                "CIDCountUnspecified",
                "CIDCountProbe",
            ]

            assay_data = {}
            # Extracting required properties from the parsed XML
            for prop in properties:
                assay_data[prop] = (
                    data_dict.get("AssaySummaries", {})
                    .get("AssaySummary", {})
                    .get(prop, None)
                )
            return assay_data
        except Exception as e:
            logging.error(f"Error parsing XML for AID {aid}: {e}")
            return None

    def extract_assay_properties(self, main_data):
        """
        Extracts detailed properties of assays from PubChem for each unique assay
        ID found in the input data file.

        This method processes an input CSV file containing assay IDs (AID) and
        performs concurrent HTTP requests to fetch detailed assay properties
        from the PubChem database. The retrieved details include assay type,
        activity name, source name, source ID, name, and description. These
        properties are compiled into a new DataFrame, which is then
        saved to a CSV file for further analysis or use.

        The method employs a ThreadPoolExecutor to manage concurrent requests
        efficiently, improving the performance when dealing with a large number
        of assay IDs. Errors encountered during data fetching are logged, and the
        process continues with the next assay ID, ensuring the method's robustness.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which
            generated after running `extractor.run()`.

        Returns:
            pd.DataFrame: A DataFrame containing the fetched assay properties,
            including columns for AID, Assay Type, Activity Name, SourceName,
            SourceID, Name, and Description. This DataFrame is saved to
            'Data/Nodes/Assay_Properties.csv' in the current working directory.

        Raises:
            ValueError: If the input CSV file is empty or does not contain the 'AID' column.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.run()
            >>> assay_properties_df = extractor.extract_assay_properties('Data/AllDataConnected.csv')
            >>> print(assay_properties_df.head())

            This example reads assay IDs from 'Data/AllDataConnected.csv',
            queries PubChem for their detailed properties, and compiles the
            results into a DataFrame, which is also saved to 'Data/Nodes/Assay_Properties.csv'.

        Note:
            This method requires network access to the PubChem API and assumes
            the availability of a valid 'AID' column in the input CSV file.
            Ensure the input file path is correct and accessible to avoid errors during processing.
        """

        df = pd.read_csv(main_data,  low_memory=False)

        # Check if the DataFrame is valid
        if df.empty or "AID" not in df.columns:
            logging.error("DataFrame is empty or does not contain 'AID' column.")
            return pd.DataFrame()

        unique_aids = df["AID"].unique().tolist()  # Extracting unique assay IDs

        columns = [
            "AID",
            "Assay Type",
            "Activity Name",
            "SourceName",
            "SourceID",
            "Name",
            "Description",
            "Protocol",
            "Comment",
            "Method",
            "Target",
            "CIDCountAll",
            "CIDCountActive",
            "CIDCountInactive",
            "CIDCountInconclusive",
            "CIDCountUnspecified",
            "CIDCountProbe",
        ]
        assay_df = pd.DataFrame(
            columns=columns
        )  # Initializing a DataFrame to store assay properties

        # Using ThreadPoolExecutor for concurrent fetching of assay details
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self._CONCURRENT_REQUEST_LIMIT
        ) as executor:
            future_to_aid = {
                executor.submit(self._fetch_assay_details, aid): aid
                for aid in unique_aids
            }

            # Iterating over completed futures
            for future in concurrent.futures.as_completed(future_to_aid):
                aid = future_to_aid[future]
                try:
                    assay_data = future.result()  # Fetching the result from the future
                    if assay_data:
                        # Preparing a new row with the fetched data
                        new_row = {
                            "AID": aid,
                            "Assay Type": df.loc[df["AID"] == aid, "Assay Type"].iloc[
                                0
                            ],
                            "Activity Name": df.loc[
                                df["AID"] == aid, "Activity Name"
                            ].iloc[0],
                            **assay_data,
                        }
                        # Adding the new row to the DataFrame
                        assay_df = pd.concat(
                            [assay_df, pd.DataFrame([new_row])], ignore_index=True
                        )
                except Exception as exc:
                    # Logging any errors encountered during the fetch
                    logging.error(f"Error occurred while processing AID {aid}: {exc}")
        retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"
        assay_df["SourceDatabase"] = "PubChem"
        assay_df["SourceEndpoint"] = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/summary/XML"
        assay_df["RetrievedAt"] = retrieval_time
        # Saving the updated DataFrame to a CSV file
        assay_df.to_csv("Data/Nodes/Assay_Properties.csv", sep=",", index=False)
        return assay_df

    def extract_protein_properties(self, main_data):
        """
        Extracts and compiles protein properties from the NCBI protein database
        based on accession numbers.

        Given a CSV file specified by `main_data`, this method reads protein
        accession numbers and performs web scraping on the NCBI protein database
        pages to extract protein titles. The method constructs a URL for
        each accession number, sends a request to retrieve the page content,
        and parses the HTML to find the protein title. The extracted titles,
        along with their corresponding accession numbers and URLs, are
        compiled into a DataFrame. This DataFrame is saved to a CSV file,
        providing a structured summary of protein properties for further analysis or use.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which
            generated after running `extractor.run()`.

        Returns:
            pd.DataFrame: A DataFrame with columns 'RefSeq Accession', 'URL',
            and 'Description', where 'Description' contains the title of the
            protein extracted from its NCBI page. This DataFrame is saved to
            'Data/Nodes/Protein_Properties.csv' in the current working directory.

        Raises:
            Exception: If there's an issue reading the initial CSV file or
            querying the NCBI database, details of the exception are logged.
            The method continues processing the next accession number,
            ensuring robustness against individual failures.

        Example:
            Assuming 'protein_data.csv' contains a column 'Target Accession'
            with accession numbers:

            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.run() # you need to run this only once
            >>> protein_properties_df = extractor.extract_protein_properties('Data/AllDataConnected.csv')
            >>> print(protein_properties_df.head())

            This would read accession numbers from 'Data/AllDataConnected.csv',
            scrape their titles from the NCBI protein database, and compile the
            results into a DataFrame, which is also saved to
            'Data/Nodes/Protein_Properties.csv'.

        Note:
            This method requires internet access to query the NCBI protein
            database. Ensure the input file path is correct and accessible to
            avoid errors during processing. Web scraping is dependent on the
            structure of the web page; changes to the NCBI protein database
            pages may require updates to the scraping logic.
        """

        # Initialize a list to store the extracted data
        data = []

        n = self._enzyme_count
        df = pd.read_csv(main_data, low_memory=False)
        gene_ids = df["Target GeneID"].value_counts().head(n).index.tolist()
        df = df[df["Target GeneID"].isin([int(item) for item in gene_ids])]
        accessions = df["Target Accession"].unique().tolist()

        for accession in accessions:
            url = f"https://www.ncbi.nlm.nih.gov/protein/{accession}"
            try:
                response = self._make_request(url)
                soup = BeautifulSoup(response.text, "html.parser")
                title = soup.title.string if soup.title else "Title Not Found"
                data.append(
                    {"RefSeq Accession": accession, "URL": url, "Description": title}
                )
            except Exception as e:
                logging.error(f"Error fetching data for accession {accession}: {e}")
                data.append(
                    {
                        "RefSeq Accession": accession,
                        "URL": url,
                        "Description": f"Error: {e}",
                    }
                )

        protein_df = pd.DataFrame(data)

        retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"
        protein_df["SourceDatabase"] = "NCBI_Protein"
        protein_df["SourceEndpoint"] = "https://www.ncbi.nlm.nih.gov/protein/{accession}"
        protein_df["RetrievedAt"] = retrieval_time

        protein_df.to_csv("Data/Nodes/Protein_Properties.csv", sep=",", index=False)
        return protein_df

    def fetch_data_batch(self, cids: List[Union[int, float]]) -> pd.DataFrame:
        """
        Retrieve detailed compound properties for a list of CIDs from PubChem
        in a single batch request (or a few chunked requests).

        Returns a DataFrame with one row per CID. Adds a StructureImage2DURL
        column per compound.
        """
        # Clean CIDs
        clean_cids = [int(c) for c in cids if pd.notna(c)]
        clean_cids = sorted(set(clean_cids))
        if not clean_cids:
            return pd.DataFrame()

        # Build the property list once
        props = (
            "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,"
            "InChIKey,IUPACName,Title,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,"
            "Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,"
            "IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,"
            "BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,"
            "PatentCount,PatentFamilyCount,LiteratureCount,Volume3D,XStericQuadrupole3D,"
            "YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,"
            "FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,"
            "FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D,"
            "Fingerprint2D"
        )

        # You can also chunk if you have thousands of CIDs
        batch_size = 150
        all_frames = []

        for i in range(0, len(clean_cids), batch_size):
            batch = clean_cids[i : i + batch_size]
            cid_str = ",".join(str(c) for c in batch)

            url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                f"{cid_str}/property/{props}/CSV"
            )

            try:
                response = self._make_request(url)
                df_batch = pd.read_csv(StringIO(response.text), sep=",", low_memory=False)

                # Add 2D structure URL per row (per CID)
                df_batch["StructureImage2DURL"] = (
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                    + df_batch["CID"].astype(str)
                    + "/PNG"
                )

                all_frames.append(df_batch)

            except Exception as e:
                logging.error(f"Error processing PubChem batch {batch[0]}–{batch[-1]}: {e}")
                continue

        if not all_frames:
            return pd.DataFrame()

        return pd.concat(all_frames, ignore_index=True)

    def fetch_data(self, cid):
        if pd.isna(cid):
            return pd.DataFrame()
        df = self.fetch_data_batch([cid])
        return df

    def extract_compound_properties(self, main_data, start_chunk=0):
        """
        Extracts and aggregates compound properties from PubChem for a list of
        compounds associated with specific genes.

        This method processes a CSV file specified by `main_data`, which contains
        gene identifiers and their associated compound IDs (CIDs). It selects
        compounds related to the top `n` most frequently occurring genes in the
        dataset, where `n` is determined by the instance's `_enzyme_count`
        attribute. The method then fetches detailed compound properties from
        PubChem in chunks, using concurrent requests to improve efficiency and
        manage the load on the PubChem API. The fetched compound properties are
        aggregated into a single DataFrame and saved to multiple CSV files,
        one for each chunk of compound IDs processed.

        Parameters:
            main_data (str): Path to a CSV file containing main data was which
            generated after running `extractor.run()`.

        Side Effects:
            - Saves the aggregated compound properties to CSV files in the current
            working directory. The files are named
            'Data/Nodes/Compound_Properties/Chunk_{i}.csv', where `{i}` is
            the chunk index.

        Returns:
            None: This method does not return a value. Instead, it saves the
            fetched compound data directly to CSV files.

        Raises:
            Exception: Logs an error and continues processing the next CID if
            an error occurs while fetching data for a specific CID.

        Example:
            >>> extractor = NodePropertiesExtractor(['CYP2D6', 'CYP3A4'])
            >>> extractor.create_data_directories()
            >>> extractor.extract_compound_properties('Data/AllDataConnected.csv')
            This will read 'Data/AllDataConnected.csv', filter for compounds
            associated with the top n genes, fetch their properties from PubChem,
            and save the results into multiple CSV files for each chunk
            of compounds processed.

        Note:
            - Ensure that the 'main_data' CSV file exists and is accessible at
            the specified path.
            - The method automatically handles NaN values in the 'CID' column
            and excludes them from processing.
            - The `enzyme_count` attribute determines the number of top genes
            for which compound properties will be fetched.
            - Internet access is required to fetch compound data from the PubChem API.
            - The method employs a `ThreadPoolExecutor` with a configurable
            number of workers (default is len(enzyme_list)) to parallelize
            requests, which can be adjusted based on system capabilities and
            API rate limits.
        """

        n = self._enzyme_count
        df = pd.read_csv(main_data, low_memory=False)
        gene_ids = df["Target GeneID"].value_counts().head(n).index.tolist()
        df = df[df["Target GeneID"].isin([int(item) for item in gene_ids])]
        df = df.dropna(subset=["CID"])
        IDs = df["CID"].unique().tolist()

        # Define chunk size and calculate number of chunks
        chunk_size = 10000
        num_chunks = math.ceil(len(IDs) / chunk_size)

        if num_chunks >= start_chunk:
            for i in range(start_chunk, num_chunks):
                # Calculate start and end indices for each chunk
                start_index = i * chunk_size
                end_index = start_index + chunk_size

                # Extract chunk of CIDs
                chunk_cids = IDs[start_index:end_index]
                # chunk_cids = [x for x in chunk_cids if not np.isnan(x)]

                # Use ThreadPoolExecutor to parallelize requests for the chunk
                with ThreadPoolExecutor(max_workers=5) as executor:
                    future_to_cid = {
                        executor.submit(self.fetch_data, cid): cid for cid in chunk_cids
                    }
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
                retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"
                chunk_df["SourceDatabase"] = "PubChem"
                chunk_df["SourceEndpoint"] = (
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/..."
                )
                chunk_df["RetrievedAt"] = retrieval_time

                # Save the concatenated DataFrame to a CSV file for the chunk
                chunk_df.to_csv(
                    f"Data/Nodes/Compound_Properties/Chunk_{i}.csv",
                    sep=",",
                    index=False,
                )
        else:
            logging.info("No more chunks to process.")


    # -------------------------------------------------------------------------
    # Derived nodes for KG schema completeness (no extra API calls)
    # -------------------------------------------------------------------------
    def extract_experimental_context_properties(
        self,
        main_data: str = "Data/AllDataConnected.csv",
        assay_properties_processed: str | None = "Data/Nodes/Assay_Properties_Processed.csv",
        gene_properties_processed: str | None = "Data/Nodes/Gene_Properties_Processed.csv",
        out_path: str = "Data/Nodes/ExperimentalContext_Properties.csv",
    ):
        """
        Build a normalized ExperimentalContext node table required by the KG schema.

        Since PubChem assay results do not reliably expose rich experimental context
        (organism, tissue/cell line, matrix, conditions) in a structured way,
        this function derives a *best-effort* context per AssayID (AID) using:
          - AssayName / AssayDescription / AssayType (when available)
          - Target Gene taxonomy (when available)

        The output is intentionally conservative: it produces one context row per assay
        with reasonable defaults and parsed hints (e.g., "human", "microsomes", "25 uM").

        Input:
          - main_data: AllDataConnected.csv (must contain AID, Assay Name, Assay Type, Target GeneID)
          - assay_properties_processed: optional, richer processed assay table (preferred)
          - gene_properties_processed: optional, to map GeneID -> TaxonomyID/Taxonomy

        Output:
          - Data/Nodes/ExperimentalContext_Properties.csv
        """
        import pandas as pd
        import numpy as np
        import os
        import re
        import datetime
        import logging

        if not os.path.exists(main_data):
            raise FileNotFoundError(f"main_data not found: {main_data}")

        df_main = pd.read_csv(main_data, low_memory=False)

        # Load processed assay properties if available (preferred source for assay name/type/description)
        df_assay = None
        if assay_properties_processed and os.path.exists(assay_properties_processed):
            try:
                df_assay = pd.read_csv(assay_properties_processed, low_memory=False)
            except Exception as exc:  # pragma: no cover
                logging.warning("Could not read %s: %s", assay_properties_processed, exc)
                df_assay = None

        # Load gene -> taxonomy mapping if available
        gene_tax = {}
        gene_tax_name = {}
        if gene_properties_processed and os.path.exists(gene_properties_processed):
            try:
                df_gene = pd.read_csv(gene_properties_processed, low_memory=False)
                if "GeneID" in df_gene.columns and "TaxonomyID" in df_gene.columns:
                    gene_tax = dict(zip(df_gene["GeneID"].astype(str), df_gene["TaxonomyID"]))
                if "GeneID" in df_gene.columns and "Taxonomy" in df_gene.columns:
                    gene_tax_name = dict(zip(df_gene["GeneID"].astype(str), df_gene["Taxonomy"]))
            except Exception as exc:  # pragma: no cover
                logging.warning("Could not read %s: %s", gene_properties_processed, exc)

        def _norm(s: str) -> str:
            return re.sub(r"\s+", " ", str(s or "")).strip().lower()

        def _infer_taxon(assay_text: str, aid: int) -> tuple[int, str]:
            t = _norm(assay_text)

            # keyword mapping
            if "human" in t or "homo sapiens" in t:
                return 9606, "Homo sapiens (human)"
            if "mouse" in t or "murine" in t:
                return 10090, "Mus musculus (mouse)"
            if "rat" in t:
                return 10116, "Rattus norvegicus (rat)"
            if "dog" in t or "canine" in t:
                return 9615, "Canis lupus familiaris (dog)"

            # fallback: most common taxonomy among genes for this assay
            sub = df_main[df_main["AID"] == aid]
            gids = sub.get("Target GeneID")
            if gids is not None:
                gids = gids.dropna().astype("Int64").astype(str).unique().tolist()
                tax_ids = [gene_tax.get(g) for g in gids if g in gene_tax]
                tax_ids = [int(x) for x in tax_ids if pd.notna(x)]
                if tax_ids:
                    # mode
                    from collections import Counter
                    tid = Counter(tax_ids).most_common(1)[0][0]
                    # try name
                    # choose any gene with that tax id
                    name = None
                    for g in gids:
                        if gene_tax.get(g) == tid:
                            name = gene_tax_name.get(g)
                            break
                    return int(tid), str(name or "Unknown")

            # default to human (common in CYP pipelines)
            return 9606, "Homo sapiens (human)"

        def _infer_matrix(assay_text: str) -> str:
            t = _norm(assay_text)
            if "microsome" in t:
                return "microsomes"
            if "hepatocyte" in t:
                return "hepatocytes"
            if "plasma" in t:
                return "plasma"
            if "serum" in t:
                return "serum"
            if "cell" in t or "cells" in t:
                return "cell culture"
            return "in vitro (unspecified)"

        def _infer_cell_tissue(assay_text: str) -> str:
            t = _norm(assay_text)
            # cell lines (small common set; extend as needed)
            for cl in ["hepg2", "hek293", "cos-7", "v79", "caco-2", "a549", "mcf-7", "hela", "u2os", "k562"]:
                if cl.lower() in t:
                    return cl.upper() if "-" not in cl else cl
            # tissues
            if "liver" in t:
                return "liver"
            if "intest" in t:
                return "intestine"
            if "kidney" in t:
                return "kidney"
            return "unspecified"

        def _infer_assay_format(assay_text: str, matrix: str) -> str:
            t = _norm(assay_text)
            if matrix == "microsomes":
                return "microsomes"
            if "cell" in t or "hepg2" in t or "hek293" in t:
                return "cellular"
            # CYP inhibition assays are often biochemical
            if "inhibition" in t or "enzyme" in t or "cyp" in t:
                return "biochemical"
            return "unknown"

        def _extract_conditions(assay_text: str) -> str:
            t = str(assay_text or "")
            # concentration patterns e.g. "25 uM"
            m = re.findall(r"(\d+(?:\.\d+)?)\s*(?:µm|um|nm|mm)\b", t, flags=re.IGNORECASE)
            concs = []
            for mm in re.finditer(r"(\d+(?:\.\d+)?)\s*(µM|uM|nM|mM)\b", t, flags=re.IGNORECASE):
                concs.append(mm.group(0))
            concs = list(dict.fromkeys(concs))
            if concs:
                return "test_concentration=" + ";".join(concs)
            return ""

        # Build per-assay contexts
        aids = pd.Series(df_main["AID"].unique()).dropna().astype(int).tolist()

        # Prepare assay metadata lookup
        assay_meta = {}
        if df_assay is not None and "AssayID" in df_assay.columns:
            cols = [c for c in ["AssayID", "AssayName", "AssayType", "AssayDescription"] if c in df_assay.columns]
            assay_meta = df_assay[cols].drop_duplicates("AssayID").set_index("AssayID").to_dict(orient="index")

        retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"
        rows = []
        for aid in aids:
            if aid in assay_meta:
                name = assay_meta[aid].get("AssayName", "")
                atype = assay_meta[aid].get("AssayType", "")
                desc = assay_meta[aid].get("AssayDescription", "")
                assay_text = " ".join([str(name), str(atype), str(desc)])
            else:
                # fallback to main_data columns
                sub = df_main[df_main["AID"] == aid]
                name = sub["Assay Name"].dropna().iloc[0] if "Assay Name" in sub.columns and sub["Assay Name"].notna().any() else ""
                atype = sub["Assay Type"].dropna().iloc[0] if "Assay Type" in sub.columns and sub["Assay Type"].notna().any() else ""
                assay_text = " ".join([str(name), str(atype)])

            tax_id, org_name = _infer_taxon(assay_text, aid)
            matrix = _infer_matrix(assay_text)
            cell_tissue = _infer_cell_tissue(assay_text)
            aformat = _infer_assay_format(assay_text, matrix)
            cond = _extract_conditions(assay_text)

            rows.append(
                {
                    "ExperimentalContextID": f"ECX:{aid}",
                    "AssayID": int(aid),
                    "OrganismName": org_name,
                    "TaxonomyID": int(tax_id) if tax_id is not None else np.nan,
                    "Matrix": matrix,
                    "CellOrTissue": cell_tissue,
                    "AssayFormat": aformat,
                    "Conditions": cond,
                    "SourceDatabase": "Derived",
                    "SourceEndpoint": main_data,
                    "RetrievedAt": retrieval_time,
                }
            )

        df_out = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        df_out.to_csv(out_path, index=False)
        logging.info("Saved ExperimentalContext properties to %s (%d rows)", out_path, len(df_out))
        return df_out


    def extract_assay_endpoint_properties(
        self,
        main_data: str = "Data/AllDataConnected.csv",
        assay_properties_processed: str | None = "Data/Nodes/Assay_Properties_Processed.csv",
        out_path: str = "Data/Nodes/AssayEndpoint_Properties.csv",
    ):
        """
        Build a normalized AssayEndpoint node table required by the KG schema.

        This table is derived from the assay result view (AllDataConnected.csv) by:
          - grouping per assay (AID) and endpoint name (Activity Name when present)
          - falling back to a qualitative endpoint ("Activity Outcome") when no endpoint name exists
          - computing basic numeric summaries over Activity Value [uM] when present

        Output:
          - Data/Nodes/AssayEndpoint_Properties.csv

        Notes:
          - AssayEndpoint IDs are stable strings: "AEP:<AID>:<endpoint_name>"
          - This node table is intentionally lightweight; BAO/EFO ontology IDs are
            added later by NodesOntologyEnricher.enrich_assay_endpoints().
        """
        import pandas as pd
        import numpy as np
        import os
        import re
        import datetime
        import logging

        if not os.path.exists(main_data):
            raise FileNotFoundError(f"main_data not found: {main_data}")

        df = pd.read_csv(main_data, low_memory=False)

        # Basic hygiene
        if "AID" not in df.columns:
            raise ValueError("main_data must contain 'AID' column")

        # Optional lookup for assay name/type to help infer readout type
        assay_meta = {}
        if assay_properties_processed and os.path.exists(assay_properties_processed):
            try:
                df_assay = pd.read_csv(assay_properties_processed, low_memory=False)
                if "AssayID" in df_assay.columns:
                    cols = [c for c in ["AssayID", "AssayName", "AssayType", "AssayDescription"] if c in df_assay.columns]
                    assay_meta = df_assay[cols].drop_duplicates("AssayID").set_index("AssayID").to_dict(orient="index")
            except Exception as exc:  # pragma: no cover
                logging.warning("Could not read %s: %s", assay_properties_processed, exc)

        def _norm(s: str) -> str:
            return re.sub(r"\s+", " ", str(s or "")).strip().lower()

        def _infer_readout_type(endpoint_name: str, aid: int) -> str:
            e = _norm(endpoint_name)
            if e in ("ic50", "ac50", "ec50", "gi50", "kc50"):
                return "potency"
            if "inhib" in e or e in ("inh",):
                return "inhibition"
            if e in ("km", "ki", "kd"):
                return "kinetic"
            # Use assay text hints
            meta = assay_meta.get(aid, {})
            assay_text = _norm(" ".join([str(meta.get("AssayName","")), str(meta.get("AssayDescription","")), str(meta.get("AssayType",""))]))
            if "induction" in assay_text:
                return "induction"
            if "inhibition" in assay_text:
                return "inhibition"
            return "unknown"

        retrieval_time = datetime.datetime.utcnow().isoformat() + "Z"

        # Normalize columns from main_data
        df["AID"] = pd.to_numeric(df["AID"], errors="coerce").astype("Int64")
        df["Activity Name"] = df.get("Activity Name")
        df["Activity Outcome"] = df.get("Activity Outcome")
        val_col = "Activity Value [uM]" if "Activity Value [uM]" in df.columns else None

        # Create endpoint_name with fallback
        endpoint_name = df["Activity Name"].astype(str)
        endpoint_name = endpoint_name.where(df["Activity Name"].notna() & (endpoint_name.str.strip() != ""), other="Activity Outcome")
        df["_EndpointName"] = endpoint_name

        group_cols = ["AID", "_EndpointName"]
        gb = df.groupby(group_cols, dropna=True)

        rows = []
        for (aid, ep), sub in gb:
            if pd.isna(aid):
                continue
            aid_int = int(aid)
            ep_str = str(ep).strip()
            if not ep_str:
                ep_str = "Activity Outcome"

            readout_type = _infer_readout_type(ep_str, aid_int)

            # Numeric summaries
            has_numeric = False
            vmin = vmed = vmax = np.nan
            unit = ""
            if val_col is not None:
                vals = pd.to_numeric(sub[val_col], errors="coerce").dropna()
                if len(vals) > 0:
                    has_numeric = True
                    vmin = float(vals.min())
                    vmed = float(vals.median())
                    vmax = float(vals.max())
                    unit = "uM"

            # Outcome counts (optional)
            oc = {}
            if "Activity Outcome" in sub.columns:
                vc = sub["Activity Outcome"].astype(str).str.strip().replace({"nan": ""})
                for k, cnt in vc.value_counts().items():
                    if not k:
                        continue
                    oc[f"OutcomeCount_{k}"] = int(cnt)

            ep_id = f"AEP:{aid_int}:{re.sub(r'[^A-Za-z0-9]+', '_', ep_str).strip('_')}"
            row = {
                "AssayEndpointID": ep_id,
                "AssayID": aid_int,
                "EndpointName": ep_str,
                "ReadoutType": readout_type,
                "Unit": unit,
                "HasNumericValue": bool(has_numeric),
                "ValueMin_uM": vmin,
                "ValueMedian_uM": vmed,
                "ValueMax_uM": vmax,
                "SourceDatabase": "Derived",
                "SourceEndpoint": main_data,
                "RetrievedAt": retrieval_time,
            }
            row.update(oc)
            rows.append(row)

        df_out = pd.DataFrame(rows)

        # Make sure every assay gets at least one endpoint row
        all_aids = pd.Series(df["AID"].dropna().astype(int).unique())
        have_aids = pd.Series(df_out["AssayID"].dropna().astype(int).unique()) if not df_out.empty else pd.Series([], dtype=int)
        missing_aids = sorted(set(all_aids.tolist()) - set(have_aids.tolist()))
        for aid in missing_aids:
            ep_str = "Activity Outcome"
            ep_id = f"AEP:{aid}:Activity_Outcome"
            df_out = pd.concat([df_out, pd.DataFrame([{
                "AssayEndpointID": ep_id,
                "AssayID": int(aid),
                "EndpointName": ep_str,
                "ReadoutType": "qualitative",
                "Unit": "",
                "HasNumericValue": False,
                "ValueMin_uM": np.nan,
                "ValueMedian_uM": np.nan,
                "ValueMax_uM": np.nan,
                "SourceDatabase": "Derived",
                "SourceEndpoint": main_data,
                "RetrievedAt": retrieval_time,
            }])], ignore_index=True)

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        df_out.to_csv(out_path, index=False)
        logging.info("Saved AssayEndpoint properties to %s (%d rows)", out_path, len(df_out))
        return df_out
