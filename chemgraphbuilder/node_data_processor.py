import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

class NodeDataPreprocessor:
    def __init__(self, data_dir: str):
        """
        Initializes the NodeDataPreprocessor with a directory path to manage the data files.

        Args:
            data_dir (str): The directory where the node data files are stored.
        """
        self.data_dir = data_dir

    def preprocess_assays(self):
        """
        Processes the assay data by renaming columns and saving the modified data back to disk.
        This method also handles visualization of assay data distributions if necessary.
        """
        df = pd.read_csv(f'{self.data_dir}/Nodes/Assay_Properties.csv')
        df.rename(columns={"AID": "AssayID", "Assay Type": "AssayType",
                           "Activity Name": "AssayActivityName", "SourceID": "AssaySourceID",
                           "SourceName": "AssaySourceName", "Name": "AssayName",
                           "Description": "AssayDescription"}, inplace=True)
        df.to_csv(f'{self.data_dir}/Nodes/Assay_Properties_Processed.csv', index=False)

    def preprocess_proteins(self):
        """
        Processes the protein data by renaming columns and saving the processed data.
        This method simplifies access to protein data for downstream analysis.
        """
        df = pd.read_csv(f'{self.data_dir}/Nodes/Protein_Properties.csv')
        df.rename(columns={"ID": "ProteinID", "Name": "ProteinName",
                           "Description": "ProteinDescription"}, inplace=True)
        df.to_csv(f'{self.data_dir}/Nodes/Protein_Properties_Processed.csv', index=False)

    def preprocess_genes(self):
        """
        Processes gene data by renaming columns and changing data types for specific fields.
        The processed data is saved for further use in gene-related analyses.
        """
        df = pd.read_csv(f'{self.data_dir}/Nodes/Gene_Properties.csv')
        df.rename(columns={"Symbol": "GeneSymbol", "Taxonomy": "TaxonomyID",
                           "Synonyms": "GeneSynonyms"}, inplace=True)
        df['GeneID'] = df['GeneID'].astype('Int64')
        df.to_csv(f'{self.data_dir}/Nodes/Gene_Properties_Processed.csv', index=False)

    def preprocess_compounds(self):
        """
        Concatenates multiple CSV files containing compound data into a single file,
        renames columns for uniformity, and saves the consolidated data. This method
        facilitates easier management and analysis of compound data.
        """
        path = f'{self.data_dir}/Nodes/Compound_Properties'
        all_csv_files = glob.glob(path + "/*.csv")
        first_file = True
        output_file = f'{path}/Compound_Properties.csv'

        with open(output_file, 'w', newline='', encoding='utf-8') as f_out:
            for file in all_csv_files:
                with open(file, 'r', newline='', encoding='utf-8') as f_in:
                    header = f_in.readline()
                    if first_file:
                        f_out.write(header)
                        first_file = False
                    for line in f_in:
                        f_out.write(line)

        df = pd.read_csv(output_file)
        df.rename(columns={"CID": "CompoundID", "Title": "CompoundName"}, inplace=True)
        df.to_csv(f'{path}/Compound_Properties_Processed.csv', index=False)
