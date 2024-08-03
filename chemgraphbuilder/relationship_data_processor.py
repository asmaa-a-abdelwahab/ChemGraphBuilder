import os
import glob
import pandas as pd
import logging
import concurrent.futures
import ast
import numpy as np

# Set up logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class RelationshipDataProcessor:
    """
    A class to process relationship data files, filtering and augmenting the data.

    Attributes:
        path (str): The directory path where the data files are stored.
        csv_files (list): List of CSV files matching the pattern 'AID_*.csv'.
        all_data_connected (dict): A dictionary containing additional data connected to assays.
    """

    def __init__(self, path, start_chunk=0):
        """
        Initializes the RelationshipDataProcessor with the specified path and start chunk index.

        Args:
            path (str): The directory path containing the CSV files.
            start_chunk (int): The starting index for processing chunks.
        """
        self.path = path
        self.csv_files = glob.glob(os.path.join(path, "AID_*.csv"))
        self.start_chunk = start_chunk
        self.all_data_connected = {}
        self.unique_column_names = []

        # Check if the all_data_connected_dict and all_columns files exist
        all_data_connected_file = 'Data/Relationships/all_data_connected_dict.txt'
        all_columns_file = 'Data/Relationships/all_columns.txt'

        if os.path.exists(all_data_connected_file):
            self.all_data_connected = self._load_all_data_connected_from_file(all_data_connected_file)
        else:
            self.all_data_connected = self._load_all_data_connected('Data/AllDataConnected.csv')

        if os.path.exists(all_columns_file):
            self.unique_column_names = self._load_columns_from_file(all_columns_file)
        else:
            self.unique_column_names = self._get_filtered_columns()
            self._save_columns_to_file(all_columns_file, self.unique_column_names)

        # Ensure the 'activity' column is included
        if 'activity' not in self.unique_column_names:
            self.unique_column_names.append('activity')

    def _load_all_data_connected(self, file_path):
        """
        Loads additional data from a specified file and organizes it into a dictionary.

        Args:
            file_path (str): The path to the file containing additional data.

        Returns:
            dict: A dictionary with keys as tuples of (aid, cid, activity_outcome)
                  and values as dictionaries of additional information.
        """
        all_data_connected = {}
        df = pd.read_csv(file_path)
        df.columns = [col.replace(' ', '_').lower() for col in df.columns]
        df = df.dropna(subset=['aid', 'cid'], how='any')
        for _, row in df.iterrows():
            key = (int(row['aid']), int(row['cid']), row['activity_outcome'])
            all_data_connected[key] = row.to_dict()

        # Optionally save the dictionary to a file
        self._save_all_data_connected_to_file(all_data_connected)

        return all_data_connected

    def _save_all_data_connected_to_file(self, all_data_connected):
        """
        Saves the all_data_connected dictionary to a file.

        Args:
            all_data_connected (dict): The dictionary to save.
        """
        with open("Data/Relationships/all_data_connected_dict.txt", "w") as file:
            for key, value in all_data_connected.items():
                file.write(f"{key}: {value}\n")

    def _load_all_data_connected_from_file(self, file_path):
        """
        Loads the all_data_connected dictionary from a file.

        Args:
            file_path (str): The path to the file containing the dictionary.

        Returns:
            dict: The loaded dictionary.
        """
        all_data_connected = {}
        with open(file_path, "r") as file:
            for line in file:
                key, value = line.strip().split(": ", 1)
                key = ast.literal_eval(key)
                # Replace 'nan' placeholder with np.nan
                value = value.replace('nan', '"__nan__"')
                value_dict = ast.literal_eval(value)
                # Convert '__nan__' back to np.nan
                for k, v in value_dict.items():
                    if v == "__nan__":
                        value_dict[k] = np.nan
                all_data_connected[key] = value_dict
        return all_data_connected

    def _get_filtered_columns(self):
        """
        Extracts unique column names from the CSV files and additional data.

        Returns:
            list: A list of unique column names.
        """
        all_columns = set()

        # Extract additional columns from the all_data_connected dictionary
        additional_columns = set()
        for value in self.all_data_connected.values():
            additional_columns.update(value.keys())

        def read_columns(file):
            try:
                # Read only column names from the CSV file
                df = pd.read_csv(file, nrows=0)
                return set([col.replace(' ', '_').lower() for col in df.columns])
            except Exception as e:
                logging.error(f"Error reading {file}: {e}")
                return set()

        # Use ThreadPoolExecutor for concurrent reading of columns from multiple files
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            results = list(executor.map(read_columns, self.csv_files))

        for columns in results:
            all_columns.update(columns)

        all_columns.update(additional_columns)

        return list(all_columns)

    def _save_columns_to_file(self, file_path, columns):
        """
        Saves the list of columns to a file.

        Args:
            file_path (str): The path to the file.
            columns (list): The list of columns to save.
        """
        with open(file_path, "w") as file:
            for item in columns:
                file.write(f"{item}\n")

    def _load_columns_from_file(self, file_path):
        """
        Loads the list of columns from a file.

        Args:
            file_path (str): The path to the file.

        Returns:
            list: The loaded list of columns.
        """
        with open(file_path, "r") as file:
            columns = [line.strip() for line in file]
        return columns

    def _add_all_data_connected_info(self, row):
        """
        Adds additional information from all_data_connected to a row.

        Args:
            row (pd.Series): A row from a DataFrame.

        Returns:
            pd.Series: The updated row with additional data if available.
        """
        key = (int(row['aid']), int(row['cid']), row['activity_outcome'])
        if key in self.all_data_connected:
            additional_info = self.all_data_connected[key]
            for col, val in additional_info.items():
                row[col] = val
        else:
            logging.warning(f"Key {key} not found in all_data_connected.")
        return row

    def process_files(self):
        """
        Processes the CSV files by filtering, cleaning, and augmenting data.

        The processed data is saved to output files.
        """
        self._filter_and_clean_data()
        logging.info("Data filtered, cleaned, and combined successfully.")

    def _filter_and_clean_data(self):
        """
        Filters and cleans data from CSV files, then saves to output files in chunks.
        """
        base_output_file = 'Data/Relationships/Assay_Compound_Relationship_Processed/Assay_Compound_Relationship'
        base_compound_gene_file = 'Data/Relationships/Compound_Gene_Relationship/Compound_Gene_Relationship'
        
        # Process files in batches of 100 
        batch_size = 100
        total_files = len(self.csv_files)

        for batch_index in range(0, total_files, batch_size):
            if batch_index >= self.start_chunk * batch_size:
                batch_files = self.csv_files[batch_index:batch_index + batch_size]
                batch_output_file = f"{base_output_file}_batch_{batch_index // batch_size + 1}.csv"
                batch_compound_gene_file = f"{base_compound_gene_file}_batch_{batch_index // batch_size + 1}.csv"

                # Initialize output files with headers for each batch
                pd.DataFrame(columns=self.unique_column_names).to_csv(batch_output_file, index=False)
                pd.DataFrame(columns=['cid', 'target_geneid', 'activity', 'aid']).to_csv(batch_compound_gene_file, index=False)

                for file in batch_files:
                    logging.info(f"Processing file {file}")
                    self._process_file(file, self.unique_column_names, batch_output_file, batch_compound_gene_file)

                logging.info(f"Processed batch {batch_index // batch_size + 1} of {total_files // batch_size + 1}")

    def _process_file(self, file, unique_column_names, output_file, compound_gene_file):
        """
        Processes a single CSV file, applying filtering, cleaning, and adding data.

        Args:
            file (str): The file path to the CSV file.
            unique_column_names (list): The list of unique column names to use.
            output_file (str): The path to the output file for combined data.
            compound_gene_file (str): The path to the output file for compound-gene relationships.
        """
        try:
            df = pd.read_csv(file, dtype={'ASSAYDATA_COMMENT': 'object'})
            df.columns = [col.replace(' ', '_').lower() for col in df.columns]
            df = df.dropna(subset=['cid'], how='any')

            phenotype_cols = [col for col in df.columns if col.startswith('phenotype')]

            if isinstance(df, pd.Series):
                df = df.to_frame().T  # Convert to DataFrame if a Series is encountered
            if df.columns.duplicated().any():
                df = df.loc[:, ~df.columns.duplicated()]
                logging.info("Duplicated columns removed from partition.")
            df = df.reindex(columns=unique_column_names, fill_value=pd.NA)
            df = df.dropna(subset=['aid', 'cid'], how='any')

            if not df.empty:
                df['measured_activity'] = df[phenotype_cols].apply(lambda row: row.mode()[0] if not row.mode().empty else None, axis=1)

                df = df.apply(self._add_all_data_connected_info, axis=1)

                if any(col in df.columns for col in phenotype_cols) and df['activity_outcome'].notna().all():
                    df = df.groupby(['activity_outcome', 'assay_name']).apply(self.propagate_phenotype).reset_index(drop=True)

                if 'target_geneid' not in df.columns:
                    df['target_geneid'] = pd.NA

                if 'sid' in df.columns:
                    df['activity_url'] = df.apply(lambda row: f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{row['aid']}#sid={row['sid']}", axis=1)
                else:
                    df['activity_url'] = pd.NA

                # Drop rows where both aid and cid are 1
                df = df[(df['aid'] != 1) | (df['cid'] != 1)]

                df = self._determine_labels_and_activity(df)

                logging.info(f"Processed file {file} with {len(df)} rows.")
                if not df.empty:
                    # Write the processed data to the output files
                    df.to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)
                    df[['cid', 'target_geneid', 'activity', 'aid']].to_csv(compound_gene_file, mode='a', header=not os.path.exists(compound_gene_file), index=False)
            else:
                logging.info(f"No data to process in file {file} after filtering.")
        except Exception as e:
            logging.error(f"Error processing file {file}: {e}")

    @staticmethod
    def most_frequent(row):
        """
        Finds the most frequent value in a row, excluding NaN values.

        Args:
            row (pd.Series): A row from a DataFrame.

        Returns:
            str: The most frequent value in the row.
        """
        values = row.dropna()
        string_values = values[values.apply(lambda x: isinstance(x, str))]
        return string_values.mode()[0] if not string_values.empty else None

    @staticmethod
    def propagate_phenotype(group):
        """
        Propagates the phenotype information within a group.

        Args:
            group (pd.DataFrame): A DataFrame group.

        Returns:
            pd.DataFrame: The updated group with propagated phenotype information.
        """
        phenotype_value = group['phenotype'].dropna().unique()
        if len(phenotype_value) > 0:
            group['phenotype'] = phenotype_value[0]
        return group

    def _determine_labels_and_activity(self, merged_df):
        """
        Determines the activity labels for the data based on predefined keywords.

        Args:
            merged_df (pd.DataFrame): The DataFrame containing merged data.

        Returns:
            pd.DataFrame: The DataFrame with determined activity labels.
        """
        inhibitor_keywords = [
            'inhibition', 'reversible inhibition', 'time dependent inhibition',
            'inhibitory activity', 'time-dependent inhibition', 'time dependent irreversible inhibition',
            'inhibitory concentration', 'inhibitory effect', 'inhibitory potency',
            'concentration required to inhibit', 'competitive inhibition', 'cyp inhibition',
            'irreversible inhibition', 'mechanism based inhibition', 'mixed inhibition',
            'mixed type inhibition', 'inhibitory constant', 'antagonistic activity', 'selectivity',
            's1p4 agonists', 'small molecule antagonists', 'displacement', 'mediated midazolam 1-hydroxylation',
            'time/nadph-dependent inhibition', 'reversal inhibition', 'mechanism-based inhibition',
            'mechanism based time dependent inhibition', 'reversible competitive inhibition',
            'predictive competitive inhibition','noncompetitive inhibition', 'in vitro inhibitory',
            'in vitro inhibition', 'inhibition of', 'direct inhibition','enzyme inhibition', 'dndi',
            'inhibition assay'
        ]

        ligand_keywords = [
            'binding affinity', 'spectral binding', 'interaction with', 'bind',
            'covalent binding affinity', 'apparent binding affinity'
        ]

        inhibitor_substrate_keywords = [
            'inhibitors and substrates'
        ]

        inhibitor_activator_modulator_keywords = [
            'apoprotein formation', 'panel assay', 'eurofins-panlabs enzyme assay'
        ]

        substrate_keywords = [
            'drug metabolism', 'prodrug', 'metabolic', 'oxidation', 'substrate activity',
            'michaelis-menten', 'metabolic stability', 'bioactivation', 'drug level',
            'enzyme-mediated drug depletion', 'enzyme-mediated compound formation',
            'phenotyping', 'activity of human recombinant cyp', 'activity of recombinant cyp',
            'activity at cyp', 'enzyme-mediated drug metabolism'
        ]

        inactivator_keywords = [
            'inactivator', 'inactivation of', 'mechanism based inactivation of', 'inactivators',
            'metabolism dependent inactivation'
        ]

        activator_keywords = [
            'assay for activators', 'activation of', 'activators of'
        ]

        inducer_keywords = [
            'induction of', 'inducer', 'inducers', 'time-dependant induction'
        ]

        all_keywords = (inhibitor_keywords + ligand_keywords + inhibitor_substrate_keywords +
                        inhibitor_activator_modulator_keywords + substrate_keywords +
                        inactivator_keywords + activator_keywords + inducer_keywords)

        keyword_to_label = {
            **{keyword: 'Inhibitor' for keyword in inhibitor_keywords},
            **{keyword: 'Inhibitor/Substrate' for keyword in inhibitor_substrate_keywords},
            **{keyword: 'Inhibitor/Inducer/Modulator' for keyword in inhibitor_activator_modulator_keywords},
            **{keyword: 'Substrate' for keyword in substrate_keywords},
            **{keyword: 'Inactivator' for keyword in inactivator_keywords},
            **{keyword: 'Activator' for keyword in activator_keywords},
            **{keyword: 'Inducer' for keyword in inducer_keywords},
            **{keyword: 'Ligand' for keyword in ligand_keywords},
        }

        def determine_active_label(assay_name):
            # Determine the appropriate label based on the first keyword found in the assay name
            assay_name_lower = assay_name.lower()
            first_keyword = None
            first_position = len(assay_name_lower)

            for keyword in all_keywords:
                position = assay_name_lower.find(keyword)
                if 0 <= position < first_position:
                    first_keyword = keyword
                    first_position = position

            if first_keyword:
                return keyword_to_label[first_keyword]
            return 'Inhibitor/Inducer/Modulator'

        merged_df['activity'] = None

        # Assign the 'Inactive' label where the activity outcome is inactive
        inactive_mask = merged_df['activity_outcome'] == 'Inactive'
        merged_df.loc[inactive_mask, 'activity'] = 'Inactive'

        # Assign labels based on assay name keywords for active outcomes
        active_mask = merged_df['activity_outcome'] == 'Active'
        if active_mask.any():
            merged_df.loc[active_mask, 'activity'] = merged_df.loc[active_mask, 'assay_name'].apply(determine_active_label)
            merged_df.loc[active_mask & merged_df['activity_name'].isin(['Km', 'Drug metabolism']), 'activity'] = 'Substrate'

            # Define the patterns with non-capturing groups
            substrate_pattern = r'(?:activity of.*oxidation)|(?:activity at cyp.*phenotyping)|(?:activity at human recombinant cyp.*formation)|(?:activity at recombinant cyp.*formation)'
            ActIndMod_pattern = r'(?:effect on cyp)|(?:effect on human recombinant cyp)|(?:effect on recombinant cyp)|(?:effect on human cyp)'
            inducer_pattern = r'(?:effect on cyp.*induction)|(?:induction of.*)'

            merged_df.loc[active_mask & merged_df['assay_name'].str.contains(substrate_pattern, case=False, regex=True), 'activity'] = 'Substrate'
            merged_df.loc[active_mask & merged_df['assay_name'].str.contains(ActIndMod_pattern, case=False, regex=True), 'activity'] = 'Inhibitor/Inducer/Modulator'
            merged_df.loc[active_mask & merged_df['assay_name'].str.contains(inducer_pattern, case=False, regex=True), 'activity'] = 'Inducer'

            merged_df.loc[active_mask & merged_df['activity_direction'].str.contains('decreasing', case=False), 'activity'] = 'Inhibitor'
            merged_df.loc[active_mask & merged_df['activity_direction'].str.contains('increasing', case=False), 'activity'] = 'Activator'
            merged_df.loc[active_mask & (merged_df['aid'] == 1215398), 'activity'] = 'Inactivator'

        return merged_df
