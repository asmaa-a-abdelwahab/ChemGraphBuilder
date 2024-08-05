"""
This module defines the `RelationshipPropertiesExtractor` class, which is responsible for extracting and analyzing
relationship properties among compounds, genes, and assays from the PubChem database.

The class facilitates the retrieval of complex relational data between chemical entities, enabling detailed analysis
of biochemical interactions and properties. The extracted data is ideal for constructing knowledge graphs, supporting
drug discovery, and understanding genetic influences on compound behavior.

Classes:
    - RelationshipPropertiesExtractor: A class to extract and analyze relationship properties from PubChem.

Usage Example:
    >>> extractor = RelationshipPropertiesExtractor()
    >>> extractor.assay_compound_relationship("Data/AllDataCollected.csv")
    This example fetches assay-compound relationship data for specified assays and saves the data to CSV files.

Note:
    Ensure network access to the PubChem API for data retrieval.
"""

import io
import os
import time
import timeit
import json
import logging
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import xml.etree.ElementTree as ET
import requests
from urllib.parse import quote, urlencode
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class RelationshipPropertiesExtractor:
    """
    Extracts and analyzes relationship properties among compounds, genes, and
    assays from the PubChem database.

    This class facilitates the retrieval of complex relational data between
    chemical entities, enabling detailed analysis of biochemical interactions
    and properties. The extracted data is ideal for constructing knowledge
    graphs, supporting drug discovery, and understanding genetic influences
    on compound behavior.

    Methods within the class are tailored to query specific relationship types
    from PubChem, including compound-assay relationships, compound co-occurrences,
    and compound transformations influenced by genes. Data fetched from PubChem
    is processed and saved in structured formats (CSV files), ready for further
    analysis or database integration.

    Attributes:
        session (requests.Session): Session object to persist certain parameters
        across requests.

    Usage:
        >>> extractor = RelationshipPropertiesExtractor()
        >>> extractor.assay_compound_relationship("Data/AllDataCollected.csv")
        This example fetches assay-compound relationship data for specified
        assays and saves the data to CSV files.
    """

    def __init__(self):
        """Initializes a RelationshipPropertiesExtractor with a Requests session
         for efficient network calls."""
        self.session = requests.Session()


    def _send_request(self, url, max_retries=5, initial_wait=1):
        for attempt in range(max_retries):
            try:
                response = self.session.get(url, timeout=30)
                response.raise_for_status()
                return response
            except requests.HTTPError as e:
                if response.status_code == 503:
                    wait = initial_wait * (2 ** attempt)
                    print(f"Server busy or under maintenance. Retrying in {wait} seconds...")
                    time.sleep(wait)
                else:
                    print(f"HTTP Error: {e}")
                    break  # Break the loop for non-503 HTTP errors
            except requests.RequestException as e:
                print(f"Request Exception: {e}")
                wait = initial_wait * (2 ** attempt)
                print(f"Network error. Retrying in {wait} seconds...")
                time.sleep(wait)
        return None  # Return None to indicate failure after all retries


    def fetch_data_for_aid(self, aid, columns_to_remove):
        """
        Fetches and processes assay data for a specified Assay ID (AID) from the
        PubChem database, preparing it for analysis or further processing.

        This method queries the PubChem database for assay data associated with
        a given AID. It constructs the query URL, sends the request using a
        previously established session, and processes the response. The response
        is expected to be in CSV format, which this method reads into a pandas
        DataFrame. Specific columns can be removed from this DataFrame based on
        the requirements for analysis. This allows for the customization of
        the fetched data, making it easier to work with specific datasets.

        If the request is successful and the data is fetched without issues,
        it undergoes initial processing to remove unwanted columns as specified
        by the 'columns_to_remove' parameter. In case of an error during the
        data fetching or processing (e.g., issues with parsing the CSV data),
        appropriate error messages are logged, and an empty DataFrame is
        returned as a fallback.

        Parameters:
            aid (int): The assay ID for which data is to be fetched. This ID is
            used to construct the query URL to the PubChem database.
            columns_to_remove (list of str): A list of column names that should
            be removed from the fetched DataFrame. This allows for the exclusion
            of data that might not be relevant to the subsequent analysis or
            processing steps.

        Returns:
            pandas.DataFrame: A DataFrame containing the processed data
            associated with the given AID. The DataFrame will exclude columns
            listed in 'columns_to_remove'. If the data fetching fails or if
            an error occurs during processing, an empty DataFrame is returned.

        Raises:
            requests.RequestException: If an error occurs during the HTTP request
            to the PubChem API. This includes scenarios such as timeout issues,
            non-200 status codes, or network-related errors. The exception is
            handled internally with logging, but it's important to be aware of
            its possibility.
            pd.errors.ParserError: If an error occurs while parsing the CSV
            response from PubChem into a DataFrame. This could happen due to
            malformed data or unexpected changes in the response format.
            Like with RequestException, this error is logged and results in
            the return of an empty DataFrame.

        Example:
            >>> extractor = RelationshipPropertiesExtractor()
            >>> processed_data_df = extractor.fetch_data_for_aid(12345, ['UnwantedColumn1', 'UnwantedColumn2'])
            >>> print(processed_data_df.head())
            This example demonstrates how to fetch and process assay data for
            the assay with ID 12345, removing 'UnwantedColumn1' and
            'UnwantedColumn2' from the resulting DataFrame. The first few rows
            of the processed DataFrame are printed as an output.

        Note:
            - This method is part of a class that requires a valid session with
            the PubChem API. Ensure that the class is properly initialized and that
            the session is active.
            - The removal of columns is an optional step and can be customized
            based on the analysis needs. If no columns need to be removed, pass an
            empty list as 'columns_to_remove'.
        """
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?"
            "query=download&record_type=datatable&actvty="
            f"all&response_type=display&aid={aid}"
        )

        response = self._send_request(url)
        if response and response.status_code == 200:
            try:
                compound_df = pd.read_csv(StringIO(response.text), sep=',')
                # Drop specified columns and process column names in-place for memory efficiency
                columns_to_remove_set = set(columns_to_remove)
                existing_columns_set = set(compound_df.columns)
                columns_to_actually_remove = list(columns_to_remove_set & existing_columns_set)
                compound_df.drop(columns=columns_to_actually_remove,
                                 errors='ignore', inplace=True)
                compound_df.rename(columns=lambda x: x.replace('PUBCHEM_', '') if x.startswith('PUBCHEM_') else x, inplace=True)

                # compound_df.drop(columns=[col for col in columns_to_remove if col in compound_df.columns], errors='ignore', inplace=True)
                # compound_df.columns = [col.replace('PUBCHEM_', '') if col.startswith('PUBCHEM_') else col for col in compound_df.columns]
                compound_df['AID'] = aid
                return compound_df
            except pd.errors.ParserError as e:
                logging.error(f"CSV parsing failed for AID {aid}: {e}")
        else:
            logging.error(f"Failed to fetch data for AID {aid}. Status code: {response.status_code if response else 'No Response'}")
        return pd.DataFrame()  # Return an empty DataFrame in case of failure



    def _process_dataframe(self, df, aid, columns_to_remove):
        """
        Processes the DataFrame by removing specified columns and renaming others.

        Parameters:
            df (pandas.DataFrame): The DataFrame to be processed.
            aid (int): The assay ID associated with the DataFrame.
            columns_to_remove (list of str): Columns to be removed from the DataFrame.
        """
        # Drop unnecessary columns efficiently
        columns_to_remove_set = set(columns_to_remove)
        df = df.drop(columns=list(columns_to_remove_set.intersection(df.columns)), errors='ignore')

        # Efficiently rename columns that start with 'PUBCHEM_'
        df.columns = [col.replace('PUBCHEM_', '') if col.startswith('PUBCHEM_') else col for col in df.columns]
        df['AID'] = aid


    def assay_compound_relationship(self, assays_data, start_chunk=0):
        """
        Processes and stores relationships between assays and compounds based
        on assay data from PubChem.

        Parameters:
            assays_data (str): Path to a CSV file containing assay IDs (AIDs).
            start_chunk (int): The starting index for processing chunks.
        """
        for chunk_idx, chunk in enumerate(pd.read_csv(assays_data, chunksize=100)):
            if chunk_idx >= start_chunk:
                columns_to_remove = ['PUBCHEM_RESULT_TAG', 'PUBCHEM_SID', 'PUBCHEM_EXT_DATASOURCE_SMILES']
                output_dir = 'Data/Relationships/Assay_Compound_Relationship'

                for aid in chunk['AID']:
                    if not os.path.exists(f'{output_dir}/AID_{aid}.csv'):
                        df = self.fetch_data_for_aid(aid, columns_to_remove)
                        if not df.empty:
                            if not os.path.exists(output_dir):
                                os.makedirs(output_dir)
                            df.to_csv(f'{output_dir}/AID_{aid}.csv', index=False)
                logging.info(f"Processed chunk {chunk_idx} for assay-compound relationships.")
            else:
                logging.info(f"No More Chunck to Process.")



    def _write_to_csv(self, df, filename):
        """
        Writes a DataFrame to a CSV file.
        """
        df.to_csv(filename, index=False)


    def assay_gene_relationship(self, main_data):
        """
        Extracts and saves relationships between assays and proteins from the
        specified dataset.

        This method processes assay data to identify relationships between
        assays and their target proteins. It selects relevant columns from the
        input data, removes duplicates to ensure unique relationships, and saves
        the cleaned data to a CSV file for further analysis or integration into
        knowledge graphs.

        Parameters:
            main_data (str): Path to the CSV file containing the main data. The
            file should include columns for 'AID' (Assay ID), 'Target GeneID',
            and 'Activity Name'.

        Returns:
            pandas.DataFrame: A DataFrame containing the unique relationships
            between assays and proteins, including the assay ID, target gene ID,
            and activity name.

        Side Effects:
            - Writes a CSV file to 'Data/Relationships/Assay_Gene_Relationship.csv',
            containing the processed relationships data.
        """
        df = pd.read_csv(main_data)
        columns_to_select = ['AID', 'Target GeneID', 'Activity Name']
        df = df[columns_to_select]
        df = df.drop_duplicates(keep='first', ignore_index=True)
        df.to_csv(f'Data/Relationships/Assay_Gene_Relationship.csv', index=False)
        return df


    def gene_protein_relationship(self, main_data):
        """
        Extracts and saves relationships between genes and proteins based on
        the provided dataset.

        This method selects relevant columns to highlight the relationships
        between genes and their corresponding proteins.
        It removes duplicate entries to ensure that each relationship is
        represented uniquely and saves the resultant data to
        a CSV file. This facilitates easy integration of genetic data into
        knowledge bases or further analysis.

        Parameters:
            main_data (str): Path to the CSV file containing gene and protein data.
            Expected columns include 'Target GeneID' and 'Target Accession'.

        Returns:
            pandas.DataFrame: A DataFrame of unique gene-protein relationships,
            including gene ID and protein accession numbers.

        Side Effects:
            - Writes the processed data to 'Data/Gene_Protein_Relationship.csv'
            in a structured CSV format.
        """
        df = pd.read_csv(main_data)
        columns_to_select = ['Target GeneID', 'Target Accession']
        df = df[columns_to_select]
        df = df.drop_duplicates(keep='first', ignore_index=True)
        df.to_csv(f'Data/Relationships/Gene_Protein_Relationship.csv', index=False)
        return df


    def compound_gene_relationship(self, main_data):
        """
        Identifies and records relationships between compounds and proteins from
        the input data.

        This method focuses on extracting compound-protein interaction data,
        including activity outcomes and values. It selects
        pertinent columns, removes duplicate records, and sorts the data by
        Compound ID and Target Accession for clarity. The cleaned dataset is
        then saved to a CSV file, providing a structured view  of how compounds
        interact with various proteins, which can be critical for drug discovery
        and pharmacological research.

        Parameters:
            main_data (str): Path to the CSV file with compound and protein data.
            This file should contain columns for 'CID' (Compound ID),
            'Target Accession', 'Activity Outcome', 'Activity Name', and
            'Activity Value [uM]'.

        Returns:
            pandas.DataFrame: A DataFrame with processed compound-protein
            relationships, sorted and cleaned for direct analysis or database
            insertion.

        Side Effects:
            - Saves the processed relationships data to
            'Data/Relationships/Compound_Gene_Relationship.csv',
            facilitating easy access and integration.
        """
        df = pd.read_csv(main_data)
        columns_to_select = ['CID', 'Target GeneID', 'Target Accession',
                             'Activity Outcome', 'Activity Name',
                             'Activity Value [uM]']
        df = df[columns_to_select]
        df = df.drop_duplicates(keep='first', ignore_index=True)
        df = df.sort_values(['CID', 'Target Accession'])
        df.dropna(axis=0 , thresh=1, inplace=True) ###
        df.to_csv(f'Data/Relationships/Compound_Gene_Relationship.csv', index=False)
        return df


    def fetch_similar_cids(self, cid):
        """
        Fetches similar compound IDs (CIDs) from the PubChem database for a
        given compound ID (CID) using 2D similarity.

        This method queries the PubChem database to find compounds that are
        similar to the given CID based on 2D structural similarity.
        The similarity threshold is set to 95%, and a maximum of 100 similar
        CIDs are fetched. The response is parsed from XML format to extract
        the similar CIDs.

        Parameters:
            cid (int): The compound ID for which similar CIDs are to be fetched.

        Returns:
            tuple: A tuple containing the original CID and a list of similar
            CIDs. If an error occurs, the list of similar CIDs will be empty.

        Raises:
            Exception: Logs an error message with the original CID and the
            exception if the request to PubChem fails or if parsing the XML
            response encounters an error.

        Note:
            - The method utilizes the `requests` library for HTTP requests and
            `xml.etree.ElementTree` for XML parsing.
            - In case of a request failure or parsing error, the method logs
            the error and returns the original CID with an empty list,
            allowing the calling function to handle the exception as needed.
        """
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
               f"fastsimilarity_2d/cid/{int(cid)}/cids/XML?Threshold=95&MaxRecords=100")
        try:
            response = requests.get(url)
            response.raise_for_status()
            xml_data = response.text

            # Parse XML data
            tree = ET.parse(io.StringIO(xml_data))
            root = tree.getroot()

            # Extracting CID values
            similar_cids = [element.text for element in root.findall('{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CID')]
            return cid, similar_cids
        except Exception as e:
            logging.error(f"Error processing CID {cid}: {e}")
            return cid, []


    def process_chunk(self, chunk):
        """
        Processes a chunk of CIDs in parallel to fetch similar CIDs for each CID
        in the chunk.

        This method uses a ThreadPoolExecutor to send out concurrent requests for
        fetching similar CIDs for a list of CIDs.
        The number of worker threads is set to 5. Each CID's request is handled
        by `fetch_similar_cids` method.

        Parameters:
            chunk (list of int): A list of compound IDs (CIDs) to process in
            parallel.

        Returns:
            list of tuples: A list of tuples, each containing a CID and its
            corresponding list of similar CIDs.

        Side Effects:
            - Utilizes concurrent threads to speed up the fetching process.
            - May log errors if any occur during the fetching of similar CIDs
            for individual CIDs.
        """
        with ThreadPoolExecutor(max_workers=5) as executor:
            futures = [executor.submit(self.fetch_similar_cids, cid) for cid in chunk]
            results = [future.result() for future in as_completed(futures)]
        return results


    def compound_similarity_relationship(self, main_data, start_chunk=0):
        """
        Identifies and records the similarity relationships between compounds
        based on a list of CIDs. The similarity is detrmined by the Tanimoto
        similarity coefficient with threshold 95% to ensure highe structural
        similarity.

        This method reads a CSV file containing compound data, filters compounds
        based on specific 'Target GeneID' values,
        and fetches similar CIDs for each compound. The compounds are processed
        in chunks to manage memory usage and improve efficiency. The results are
        saved into separate CSV files for each chunk.

        Parameters:
            main_data (str): Path to the CSV file containing the main compound data.
            start_chunk (int): The starting index for processing chunks.
        Note:
            - The method filters the main data for compounds associated with
            specific 'Target GeneID' values before fetching similar CIDs,
            optimizing the process for relevant compounds only.
            - The division of CIDs into chunks and concurrent processing helps
            in managing large datasets and utilizes parallelism for faster
            execution.
        """
        df = pd.read_csv(main_data)
        df = df[df['Target GeneID'].isin([1576, 1544, 1557, 1559, 1565])]
        df = df.dropna(subset=['CID'])
        IDs = df['CID'].unique().tolist()

        chunk_size=10000
        chunks = [IDs[i:i + chunk_size] for i in range(0, len(IDs), chunk_size)]

        for i, chunk in enumerate(chunks, start=0):
            if i >= start_chunk:
                chunk_results = self.process_chunk(chunk)
                chunk_df = pd.DataFrame(chunk_results, columns=['CID', 'Similar CIDs'])
                if not os.path.exists('Data/Relationships/Compound_Similarities'):
                    os.makedirs('Data/Relationships/Compound_Similarities')
                chunk_df.to_csv(f'Data/Relationships/Compound_Similarities/Chunk_{i}.csv', index=False)
                logging.info(f"Processed chunk {i} for compound similarity relationships.")


    def _fetch_chemical_neighbor_data(self, cid):
        """
        Fetches chemical-chemical relationship data for a given CID.

        Args:
            cid (int): The compound ID for which data is to be fetched.

        Returns:
            list: List of chemical-chemical relationship data.
        """
        cpd_cpd_url = ("https://pubchem.ncbi.nlm.nih.gov/link_db/link_db_server.cgi?format=JSON&type="
                       f"ChemicalNeighbor&operation=GetAllLinks&id_1={int(cid)}&response_type=display")
        # print('url', cpd_cpd_url)
        try:
            response = self._send_request(cpd_cpd_url)
            data = response.json()
            # print(data)
            return data.get('LinkDataSet', {}).get('LinkData', [])
        except Exception as e:
            logging.error(f"Failed to fetch chemical-chemical data for CID {cid}: {e}")
            return []


    def _fetch_chemical_gene_data(self, gid):
        """
        Fetches chemical-gene relationship data for a given CID.

        Args:
            cid (int): The compound ID for which data is to be fetched.

        Returns:
            list: List of chemical-gene relationship data.
        """
        cpd_gene_url = ("https://pubchem.ncbi.nlm.nih.gov/link_db/link_db_server.cgi?format=JSON&"
                        f"type=GeneSymbolChemicalNeighbor&operation=GetAllLinks&id_1={gid}&response_type=display")
        try:
            response = self._send_request(cpd_gene_url)
            data = response.json()
            return data.get('LinkDataSet', {}).get('LinkData', [])
        except Exception as e:
            logging.error(f"Failed to fetch chemical-gene data for CID {gid}: {e}")
            return []

    
    def _fetch_chemical_gene_interaction_data(self, gid):
        """
        Fetches chemical-gene relationship data for a given CID.

        Args:
            cid (int): The compound ID for which data is to be fetched.

        Returns:
            list: List of chemical-gene relationship data.
        """
        base_url = "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi"
        query = {
            "download": "*",
            "collection": "consolidatedcompoundtarget",
            "order": ["cid,asc"],
            "start": 1,
            "limit": 10000000,
            "downloadfilename": f"pubchem_geneid_{int(gid)}_consolidatedcompoundtarget",
            "where": {
                "ands": [
                    {"geneid": f"{int(gid)}"}
                ]
            }
        }
        
        # Convert the query dictionary to a JSON string
        query_json_str = json.dumps(query)
        
        # URL encode the JSON string
        encoded_query = quote(query_json_str)
        
        # Construct the full URL
        url = f"{base_url}?infmt=json&outfmt=json&query={encoded_query}"

        try:
            response = self._send_request(url)
            data = response.json()
            return data
        except Exception as e:
            logging.error(f"Failed to fetch chemical-gene data for CID {gid}: {e}")
            return []
            

    def _write_data_to_csv(self, data, filename, filter_condition=None):
        """
        Writes given data to a CSV file, with optional filtering before saving.

        This method takes a list of dictionaries (data), converts it into a
        pandas DataFrame, and optionally filters the DataFrame based on
        specified conditions before writing the result to a CSV file. The
        filtering is performed on specified columns with their expected
        values provided in 'filter_condition'. This allows for selective
        data saving, especially useful when dealing with large datasets
        or when only a subset of data is needed for further processing
        or analysis.

        Parameters:
            data (list of dict): Data to be written to a CSV file. Each
            dictionary in the list represents a row in the DataFrame, with keys
            as column names and values as row values.
            filename (str): Path to the CSV file where the data will be saved.
            If the file exists, it will be overwritten.
            filter_condition (dict, optional): A dictionary specifying the
            columns to filter by and the values to include. Keys in the
            dictionary are column names, and values are lists of acceptable
            values for that column. Rows not meeting the filter condition are
            excluded from the final DataFrame to be saved.

        Side Effects:
            - Writes a CSV file to the given filename path. The file is overwritten
            if it already exists.
            - Logs a warning if a specified column for filtering is not found in
            the DataFrame.
        """

        df = pd.DataFrame(data)
        if filter_condition:
            for column, values in filter_condition.items():
                if column in df.columns:
                    df = df[df[column].isin(values)]
                else:
                    logging.warning(f"Column {column} not found in DataFrame.")
        if not df.empty:
            df.to_csv(filename, index=False)

        
    def compound_compound_cooccurrence(self, main_data, rate_limit=5):
        """
        Analyzes compound-compound co-occurrence relationships from the specified main data file and saves the results into structured CSV files.
        
        Args:
            main_data (str): Path to the main data file.
            rate_limit (int): The maximum number of requests per second.
        
        Returns:
            str: A message indicating the completion of data fetching and saving.
        """
        logging.info("Starting compound-compound co-occurrence data retrieval ...")
        start_time = timeit.default_timer()
        
        try:
            df = pd.read_csv(main_data)
            logging.info(f"Loaded data from {main_data}. Total rows: {len(df)}")
        except FileNotFoundError:
            logging.error(f"File not found: {main_data}")
            return "File not found."
        except pd.errors.EmptyDataError:
            logging.error(f"Empty data file: {main_data}")
            return "Empty data file."
        except Exception as e:
            logging.error(f"Error reading {main_data}: {e}")
            return "Error reading data file."
    
        compound_ids = df['CID'].dropna().unique().tolist()
        logging.info(f"Unique Compound IDs to process: {len(compound_ids)}")
    
        for compound_id in compound_ids:
            logging.info(f"Processing Compound ID {int(compound_id)}")
            try:
                data = self._fetch_chemical_neighbor_data(int(compound_id))
                filename = f"Data/Relationships/Cpd_Cpd_CoOccurrence/CID_{int(compound_id)}.csv"
                self._write_data_to_csv(data, filename)
                logging.info(f"Successfully wrote data for Compound ID {int(compound_id)} to {filename}")
            except Exception as e:
                logging.error(f"Error processing Compound ID {int(compound_id)}: {e}")
            time.sleep(1 / rate_limit)  # Ensuring we don't exceed rate limit
    
        elapsed = timeit.default_timer() - start_time
        logging.info(f"Compound-compound data fetching and saving completed in {elapsed:.2f} seconds.")
        return "Compound-compound data fetching and saving completed."


    def compound_gene_cooccurrence(self, gene_data, rate_limit=5):
        """
        Analyzes compound-gene co-occurrence relationships from the specified main data file and saves the results into structured CSV files.
        """
        logging.info("Starting compound-gene co-occurrence analysis...")
        start_time = timeit.default_timer()
        
        try:
            df = pd.read_csv(gene_data)
            logging.info(f"Loaded data from {gene_data}. Total rows: {len(df)}")
        except FileNotFoundError:
            logging.error(f"File not found: {gene_data}")
            return "File not found."
        except pd.errors.EmptyDataError:
            logging.error(f"Empty data file: {gene_data}")
            return "Empty data file."
        except Exception as e:
            logging.error(f"Error reading {gene_data}: {e}")
            return "Error reading data file."

        gene_symbols = df['GeneSymbol'].unique().tolist()
        logging.info(f"Unique Gene Symbols to process: {len(gene_symbols)}")

        for gene_symbol in gene_symbols:
            logging.info(f"Processing Gene Symbol {gene_symbol}")
            try:
                data = self._fetch_chemical_gene_data(gene_symbol)
                filename = f"Data/Relationships/Cpd_Gene_CoOccurrence/Cpd_Gene_CoOccurrence_{gene_symbol}.csv"
                self._write_data_to_csv(data, filename)
                logging.info(f"Successfully wrote data for Gene Symbol {gene_symbol} to {filename}")
            except Exception as e:
                logging.error(f"Error processing Gene Symbol {gene_symbol}: {e}")
            time.sleep(1 / rate_limit)  # Ensuring we don't exceed rate limit

        elapsed = timeit.default_timer() - start_time
        logging.info(f"Compound-gene data fetching and saving completed in {elapsed:.2f} seconds.")
        return "Compound-gene data fetching and saving completed."

    
    def compound_gene_interaction(self, gene_data, rate_limit=5):
        """
        Analyzes compound-gene co-occurrence relationships from the specified main data file and saves the results into structured CSV files.
        """
        logging.info("Starting compound-gene co-occurrence analysis...")
        start_time = timeit.default_timer()
        
        try:
            df = pd.read_csv(gene_data)
            logging.info(f"Loaded data from {gene_data}. Total rows: {len(df)}")
        except FileNotFoundError:
            logging.error(f"File not found: {gene_data}")
            return "File not found."
        except pd.errors.EmptyDataError:
            logging.error(f"Empty data file: {gene_data}")
            return "Empty data file."
        except Exception as e:
            logging.error(f"Error reading {gene_data}: {e}")
            return "Error reading data file."

        gene_symbols = df['GeneID'].unique().tolist()
        logging.info(f"Unique Gene Symbols to process: {len(gene_symbols)}")

        for gene_symbol in gene_symbols:
            logging.info(f"Processing Gene Symbol {gene_symbol}")
            try:
                data = self._fetch_chemical_gene_interaction_data(int(gene_symbol))
                filename = f"Data/Relationships/Compound_Gene_Relationship/Compound_Gene_Interaction_Outside_PubChem_{int(gene_symbol)}.csv"
                df = pd.DataFrame(data)
                if not df.empty:
                    # Reorder Columns
                    all_columns = [col for col in df.columns if col not in ('cid', 'geneid')]
                    all_columns = ['cid', 'geneid'] + all_columns
                    df = df[all_columns]
                    df.to_csv(filename, index=False)
                # self._write_data_to_csv(data, filename)
                logging.info(f"Successfully wrote data for Gene Symbol {int(gene_symbol)} to {filename}")
            except Exception as e:
                logging.error(f"Error processing Gene Symbol {int(gene_symbol)}: {e}")
            time.sleep(1 / rate_limit)  # Ensuring we don't exceed rate limit

        elapsed = timeit.default_timer() - start_time
        logging.info(f"Compound-gene data fetching and saving completed in {elapsed:.2f} seconds.")


    def compound_transformation(self, gene_properties):
        """
        Analyzes compound transformation data based on gene properties, focusing
        on metabolic transformations involving specified genes. This method
        queries the PubChem database for transformation data related
        to compounds associated with the genes identified in the provided CSV file.
    
        Parameters:
            gene_properties (str): Path to the CSV file containing gene properties
            generated by the NodePropertiesExtractor class, which should include
            'GeneID' as one of its columns. This file is used to identify genes
            of interest for which compound transformation data will be fetched.
    
        Processing Steps:
            1. Reads the provided CSV file to extract unique gene identifiers.
            2. For each gene identifier, constructs a query to fetch relevant
            compound transformation data from PubChem, focusing on metabolic
            transformations where the gene plays a role.
            3. Processes and aggregates the fetched data into a structured
            pandas DataFrame.
            4. Filters the aggregated data to retain specific columns relevant
            to compound transformations, including substrate and metabolite
            Compound IDs (CIDs), the type of metabolic conversion, gene
            identifiers, PubMed IDs, and DOIs for related publications.
            5. Saves the aggregated and filtered DataFrame to a CSV file for
            further analysis or integration into knowledge graphs or other
            data models.
    
        Returns:
            pandas.DataFrame: A DataFrame containing processed compound
            transformation data, including substrate and metabolite CIDs,
            metabolic conversion types, gene identifiers, PubMed IDs, and DOIs.
            The DataFrame structure facilitates further analysis or use in
            constructing detailed views of metabolic pathways involving the
            specified genes.
    
        Side Effects:
            - Saves the aggregated compound transformation data to
            'Data/Relationships/Compound_Transformation.csv'
            in the current working directory. This file captures the relationship
            between substrates, metabolites, and genes based on the input gene
            properties.
    
        Raises:
            FileNotFoundError: If the specified 'gene_properties' file does not
            exist or cannot be read.
            ValueError: If 'gene_properties' does not contain the required
            'GeneID' column.
    
        Example:
            >>> extractor = RelationshipPropertiesExtractor()
            >>> transformation_df = extractor.compound_transformation('Data/Nodes/gene_properties.csv')
            >>> print(transformation_df.head())
            This example processes gene properties from
            'path/to/gene_properties.csv', queries PubChem for
            compound transformation data related to the genes,
            and compiles the results into a DataFrame.
    
        Note:
            The method assumes that the input 'gene_properties' file is
            accessible and correctly formatted.
            The availability and structure of the PubChem database may affect
            the completeness and accuracy of the fetched transformation data.
            Users should verify the existence of the 'Data/Relationships'
            directory and have appropriate permissions to write files to it.
        """
        df = pd.read_csv(gene_properties)
        IDs = df['Target GeneID'].unique().tolist()
    
        transformation_dfs = []
    
        for gid in IDs:
            if not np.isnan(gid):
                gid = int(gid)
                base_url = 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query='
                query = {
                    "download": "*",
                    "collection": "chemblmetabolism",
                    "order": ["relevancescore,desc"],
                    "start": 1,
                    "limit": 10000000,
                    "downloadfilename": f"pubchem_geneid_{gid}_chemblmetabolism",
                    "where": {
                        "ands": [{"geneid": gid}]
                    }
                }
                
                # Convert the dictionary to a JSON string
                query_string = json.dumps(query)
                # URL encode the JSON string
                encoded_query = quote(query_string)
                # Construct the final URL
                url = f"{base_url}{encoded_query}"
    
                response = self._send_request(url)
                if response:
                    try:
                        # Read the CSV data
                        transformation_df = pd.read_csv(StringIO(response.text), sep=',', header=0, low_memory=False)
                        print(response.text)
    
                        # Ensure columns exist
                        transformation_df = transformation_df[['substratecid',
                                                               'metabolitecid',
                                                               'metconversion',
                                                               'geneids',
                                                               'pmids',
                                                               'dois']]
    
                        # Append the DataFrame to the list
                        transformation_dfs.append(transformation_df)
                    except pd.errors.ParserError as e:
                        logging.error(f"Error parsing CSV for gene ID {gid}: {e}\nurl:{url}")
                        continue  # Skip this gene ID and continue with others
    
        # Concatenate all DataFrames
        if transformation_dfs:
            transformation_df = pd.concat(transformation_dfs, ignore_index=True)
        else:
            transformation_df = pd.DataFrame(columns=['substratecid', 'metabolitecid', 'metconversion', 'geneids', 'pmids', 'dois'])
    
        self._write_to_csv(transformation_df, 'Data/Relationships/Compound_Transformation.csv')
    
        return transformation_df
