"""
node_data_processor.py

This module provides the NodeDataProcessor class, which is responsible for
preprocessing various types of node data (assays, proteins, genes, and compounds)
for use in chemical knowledge graph construction. The preprocessing includes
renaming columns, consolidating multiple files, and saving the processed data
in a consistent format. This step ensures uniformity and ease of access for
subsequent data analysis and integration processes.

Classes:
    NodeDataProcessor: Handles preprocessing of assay, protein, gene, and compound data.

Example Usage:
    >>> processor = NodeDataProcessor(data_dir='path/to/data')
    >>> processor.preprocess_assays()
    >>> processor.preprocess_proteins()
    >>> processor.preprocess_genes()
    >>> processor.preprocess_compounds()
"""

import glob
import pandas as pd
import logging

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class NodeDataProcessor:
    """
    NodeDataProcessor is responsible for preprocessing various types of node data
    (assays, proteins, genes, and compounds) by renaming columns, consolidating
    multiple files, and saving the processed data. This preprocessing step is
    crucial for ensuring uniformity and ease of access in subsequent analysis
    and integration processes.

    Attributes:
        data_dir (str): The directory where the node data files are stored.

    Methods:
        preprocess_assays(): Processes and renames columns in assay data.
        preprocess_proteins(): Processes and renames columns in protein data.
        preprocess_genes(): Processes and renames columns in gene data.
        preprocess_compounds(): Consolidates and renames columns in compound data.
    """

    def __init__(self, data_dir: str):
        """
        Initializes the NodeDataProcessor with a directory path to manage the data files.

        Args:
            data_dir (str): The directory where the node data files are stored.
        """
        self.data_dir = data_dir

    def preprocess_assays(self):
        """
        Processes the assay data by renaming columns and saving the modified data back to disk.
        This method also handles visualization of assay data distributions if necessary.
        """
        df = pd.read_csv(f"{self.data_dir}/Nodes/Assay_Properties.csv")
        df.rename(
            columns={
                "AID": "AssayID",
                "Assay Type": "AssayType",
                "Activity Name": "AssayActivityName",
                "SourceID": "AssaySourceID",
                "SourceName": "AssaySourceName",
                "Name": "AssayName",
                "Description": "AssayDescription",
            },
            inplace=True,
        )
        df.to_csv(f"{self.data_dir}/Nodes/Assay_Properties_Processed.csv", index=False)

    def preprocess_proteins(self):
        """
        Processes the protein data by renaming columns and saving the processed data.
        This method simplifies access to protein data for downstream analysis.
        """
        df = pd.read_csv(f"{self.data_dir}/Nodes/Protein_Properties.csv")
        df.rename(
            columns={
                "RefSeq Accession": "ProteinRefSeqAccession",
                "Description": "ProteinDescription",
            },
            inplace=True,
        )
        df.to_csv(
            f"{self.data_dir}/Nodes/Protein_Properties_Processed.csv", index=False
        )

    def preprocess_genes(self):
        """
        Processes gene data by renaming columns and changing data types for specific fields.
        The processed data is saved for further use in gene-related analyses.
        """
        df = pd.read_csv(f"{self.data_dir}/Nodes/Gene_Properties.csv")
        df.rename(
            columns={
                "Symbol": "GeneSymbol",
                "Taxonomy ID": "TaxonomyID",
                "Synonyms": "GeneSynonyms",
            },
            inplace=True,
        )
        df["GeneID"] = df["GeneID"].astype("Int64")
        df["TaxonomyID"] = df["TaxonomyID"].astype("Int64")
        df.to_csv(f"{self.data_dir}/Nodes/Gene_Properties_Processed.csv", index=False)

    def preprocess_compounds(self):
        """
        Concatenates multiple CSV files containing compound data into a single file,
        renames columns for uniformity, and saves the consolidated data. This method
        facilitates easier management and analysis of compound data.
        """
        path = f"{self.data_dir}/Nodes/Compound_Properties"
        all_csv_files = glob.glob(path + "/*.csv")
        first_file = True
        output_file = "Data/Nodes/Compound_Properties.csv"

        with open(output_file, "w", newline="", encoding="utf-8") as f_out:
            for file in all_csv_files:
                with open(file, "r", newline="", encoding="utf-8") as f_in:
                    header = f_in.readline()
                    if first_file:
                        f_out.write(header)
                        first_file = False
                    for line in f_in:
                        f_out.write(line)

        df = pd.read_csv(output_file)
        df.rename(columns={"CID": "CompoundID", "Title": "CompoundName"}, inplace=True)
        df.to_csv(f"{output_file.replace('.csv', '_Processed')}", index=False)
