"""
relationship_data_processor.py

This module provides the RelationshipDataProcessor class, which processes
relationship data from multiple CSV files, cleans the data, merges it into a
single CSV file, and performs additional data processing.

The primary purpose of this module is to facilitate the extraction, cleaning,
and merging of assay-compound relationship data from various sources into a
cohesive dataset that can be used for further analysis and research.

Classes:
    RelationshipDataProcessor: Handles the processing and merging of
    relationship data from multiple CSV files.

Example Usage:
    >>> processor = RelationshipDataProcessor('/path/to/data')
    >>> processor.process_files()
"""

import os
import glob
import pandas as pd

class RelationshipDataProcessor:
    """
    A class to process relationship data from multiple CSV files, clean the data,
    merge it into a single CSV file, and perform additional data processing.

    Attributes:
    path (str): The path to the directory containing the CSV files.
    output_files (dict): A dictionary to store unique headers and corresponding 
    output file paths.
    csv_files (list): A list of all CSV files to be processed.
    processed_csv_files (list): A list of processed CSV files.
    """

    def __init__(self, path):
        """
        Initialize the RelationshipDataProcessor with the specified path.

        Parameters:
        path (str): The path to the directory containing the CSV files.
        """
        self.path = path
        self.output_files = {}
        self.csv_files = glob.glob(os.path.join(path, "AID_*.csv"))
        self.processed_csv_files = glob.glob(os.path.join(path,
                                                          "Assay_Compound_Relationship*.csv"))

    def process_files(self):
        """
        Process all CSV files by reading, cleaning, merging, and formatting the data.
        """
        # self._process_individual_files()
        # print("Step 1: Individual files processed and cleaned.")

        # self._filter_and_combine_data()
        # print("Step 2: Data filtered and combined successfully.")

        self._clean_and_save_final_data()
        print("Step 3: Final data processing and merging completed.")

    def _process_individual_files(self):
        """
        Process all individual CSV files: read, clean, and save them based on their headers.
        """
        for file in self.csv_files:
            with open(file, 'r', encoding='utf-8') as f:
                header = next(f).strip()
                num_columns = len(header.split(','))

                if header not in self.output_files:
                    output_filename = f"{self.path}/Assay_Compound_Relationship{len(self.output_files)}_cols-{num_columns}.csv"
                    self.output_files[header] = output_filename
                    with open(output_filename, 'w', encoding='utf-8') as f_out:
                        f_out.write(header + '\n')

            df = pd.read_csv(file)
            df = df[header.split(',')]
            df.dropna(subset=['CID'], inplace=True)
            df.dropna(axis=1, how='all', inplace=True)
            df.to_csv(self.output_files[header], mode='a', index=False, header=False)

    def _filter_and_combine_data(self):
        """
        Filter and format the merged data by selecting specific columns and 
        saving it to a final output file.
        """

        output_file = os.path.join(self.path, 'Assay_Compound_DataSet.csv')

        unique_column_names = self._get_filtered_columns()

        # Create and write the header of the output file if not exists
        file_exists = os.path.isfile(output_file)
        if file_exists:
            os.remove(output_file)
        pd.DataFrame(columns=unique_column_names).to_csv(output_file, index=False)

        # Read each file, adjust columns, and append to the CSV
        for file in self.processed_csv_files:
            try:
                # Load just the column names first
                df_headers = pd.read_csv(file, nrows=0)
                df_headers.columns = [col.lower() for col in df_headers.columns]

                # Check and handle duplicate column names
                if df_headers.columns.duplicated().any():
                    df_headers = df_headers.loc[:, ~df_headers.columns.duplicated()]
                    print(f"Duplicated columns removed from {file}")

                # Load the data in chunks using the cleaned column names
                for df in pd.read_csv(file, chunksize=10000,
                                      names=df_headers.columns, header=0):
                    # Reindex to the union of columns, filling missing ones with NaNs
                    df = df.reindex(columns=unique_column_names)
                    # Append to the output file
                    df.to_csv(output_file, mode='a', header=False, index=False)

                print(f"Successfully processed and appended data from {file}")
            except Exception as e:
                print(f"Error processing file {file}: {e}")

        print(f"Data appended successfully to {output_file}")

    def _get_filtered_columns(self):
        """
        Determine the columns to be included in the final output file based on 
        specific keywords.

        Returns:
        list: A list of filtered columns.
        """
        all_columns = set()
        for file in self.processed_csv_files:
            df = pd.read_csv(file, nrows=0)
            df.columns = [col.lower() for col in df.columns]
            all_columns.update(df.columns)

        keywords = ['id', 'outcome', 'phenotype', 'activity_url']
        filtered_columns = sorted([col for col in all_columns if any(keyword in col for keyword in keywords)])
        print(f"Filtered columns: {filtered_columns}")
        return filtered_columns

    def _get_filtered_columns(self):
        """
        Determine the columns to be included in the final output file based on specific keywords.

        Returns:
        list: A list of filtered columns.
        """
        all_columns = set()
        for file in self.processed_csv_files:
            df = pd.read_csv(file, nrows=0)
            df.columns = [col.lower() for col in df.columns]
            all_columns.update(df.columns)

        keywords = ['id', 'outcome', 'phenotype', 'activity_url','activity direction']
        filtered_columns = sorted([col for col in all_columns if any(keyword in col for keyword in keywords)])
        print(f"Filtered columns: {filtered_columns}")
        return filtered_columns

    def _clean_and_save_final_data(self):
        """
        Clean and save the processed dataset.
        """
        final_file_path = os.path.join(self.path, 'Assay_Compound_DataSet.csv')
        if os.path.exists(final_file_path):
            dataset = pd.read_csv(final_file_path)
            # print(dataset.columns)
            dataset.dropna(subset=['aid', 'cid'], how='any', axis=0, inplace=True)

            phenotype_cols = [col for col in dataset.columns if col.startswith('phenotype')]
            dataset['measured_activity'] = dataset[phenotype_cols].apply(self.most_frequent, axis=1)
            dataset = dataset.drop(columns=phenotype_cols)
            dataset.rename(columns={
                'activity_outcome': 'Activity Outcome',
                'activity_url': 'Activity URL',
                'measured_activity': 'Phenotype',
                'aid': 'AID',
                'cid': 'CID'
            }, inplace=True)

            AllDataConnected = pd.read_csv('Data/AllDataConnected.csv')

            # Process the dataset
            merged_df = dataset.merge(AllDataConnected,
                                      on=['AID', 'CID', 'Activity Outcome'],
                                      how='left')
            merged_df.drop(['PubMed ID', 'RNAi'], axis=1, inplace=True)
            merged_df.dropna(axis=1, how='all', inplace=True)
            merged_df = merged_df.groupby(['Activity Outcome', 'Assay Name']).apply(self.propagate_phenotype).reset_index(drop=True)
            merged_df = merged_df[merged_df['Target GeneID'].isin([1576, 1565, 1559, 1557, 1544])]
            merged_df['AID'] = merged_df['AID'].astype(int)
            merged_df['CID'] = merged_df['CID'].astype(int)
            merged_df['Target GeneID'] = merged_df['Target GeneID'].astype(int)
            merged_df['Activity URL'] = merged_df.apply(lambda row: f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{row['AID']}#sid={row['SID']}", axis=1)
            # Additional processing for determining labels and activity
            merged_df = self._determine_labels_and_activity(merged_df)
            merged_df.to_csv('Data/Relationships/Assay_Compound_Relationship/Assay_Compound_DataSet_Processed.csv', index=False)

        else:
            print(f"File {final_file_path} does not exist. No data to process.")

    @staticmethod
    def most_frequent(row):
        """
        Find the most frequent string value in a row, excluding NaNs and 
        non-string types.

        Parameters:
        row (Series): The row of data.

        Returns:
        str or None: The most frequent string value, or None if no strings are found.
        """
        values = row.dropna()
        string_values = values[values.apply(lambda x: isinstance(x, str))]
        return string_values.mode()[0] if not string_values.empty else None

    @staticmethod
    def propagate_phenotype(group):
        """
        Propagate non-null values within each group for the 'Phenotype' column.

        Parameters:
        group (DataFrame): The group of data.

        Returns:
        DataFrame: The group with propagated 'Phenotype' values.
        """
        phenotype_value = group['Phenotype'].dropna().unique()
        if len(phenotype_value) > 0:
            group['Phenotype'] = phenotype_value[0]
        return group

    def _determine_labels_and_activity(self, merged_df):
        """
        Determine the labels and activity based on keywords in the assay names.
        """

        # Define the mapping values
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
            'inactivator', 'inactivation of', 'mechanism based inactivation of',
            'inactivators', 'metabolism dependent inactivation'
        ]

        activator_keywords = [
            'assay for activators', 'activation of', 'activators of'
        ]

        inducer_keywords = [
            'induction of', 'inducer', 'inducers', 'time-dependant induction'
        ]

        # Combine all keywords into a single list for easy checking
        all_keywords = (inhibitor_keywords + ligand_keywords + inhibitor_substrate_keywords +
                        inhibitor_activator_modulator_keywords + substrate_keywords +
                        inactivator_keywords + activator_keywords + inducer_keywords)

        # Define a dictionary to map each keyword to its corresponding label
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

        # Function to determine the appropriate label based on the first keyword appearance in the assay name
        def determine_active_label(assay_name):
            assay_name_lower = assay_name.lower()
            first_keyword = None
            first_position = len(assay_name_lower)  # Start with the maximum length of the assay name

            for keyword in all_keywords:
                position = assay_name_lower.find(keyword)
                if 0 <= position < first_position:  # Check if keyword is found and is the earliest one so far
                    first_keyword = keyword
                    first_position = position

            if first_keyword:
                return keyword_to_label[first_keyword]
            return 'Inhibitor/Inducer/Modulator'

        # Set default activity to None
        merged_df['Activity'] = None

        # Handle Inactive outcomes
        inactive_mask = merged_df['Activity Outcome'] == 'Inactive'
        merged_df.loc[inactive_mask, 'Activity'] = 'Inactive'

        # Handle Active outcomes
        active_mask = merged_df['Activity Outcome'] == 'Active'
        if active_mask.any():
            merged_df.loc[active_mask, 'Activity'] = merged_df.loc[active_mask, 'Assay Name'].apply(determine_active_label)
            merged_df.loc[active_mask & merged_df['Activity Name'].isin(['Km', 'Drug metabolism']), 'Activity'] = 'Substrate'

            # Apply the changes to the 'Activity' column based on the combined conditions
            substrate_pattern = r'(activity of.*oxidation)|(activity at cyp.*phenotyping)|(activity at human recombinant cyp.*formation)|(activity at recombinant cyp.*formation)'
            merged_df.loc[active_mask & merged_df['Assay Name'].str.contains(substrate_pattern, case=False, regex=True), 'Activity'] = 'Substrate'

            ActIndMod_pattern = r'(effect on cyp)|(effect on human recombinant cyp)|(effect on recombinant cyp)|(effect on human cyp)'
            merged_df.loc[active_mask & merged_df['Assay Name'].str.contains(ActIndMod_pattern, case=False, regex=True), 'Activity'] = 'Inhibitor/Inducer/Modulator'

            inducer_pattern = r'(effect on cyp.*induction)|(induction of.*)'
            merged_df.loc[active_mask & merged_df['Assay Name'].str.contains(inducer_pattern, case=False, regex=True), 'Activity'] = 'Inducer'

            merged_df.loc[active_mask & merged_df['activity direction'].str.contains('decreasing', case=False), 'Activity'] = 'Inhibitor'
            merged_df.loc[active_mask & merged_df['activity direction'].str.contains('increasing', case=False), 'Activity'] = 'Activator'
            merged_df.loc[active_mask & (merged_df['AID'] == 1215398), 'Activity'] = 'Inactivator'

        return merged_df
