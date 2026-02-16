"""
node_data_processor.py

This module provides the NodeDataProcessor class, which is responsible for
preprocessing various types of node data (assays, proteins, genes, compounds, experimental contexts, and assay endpoints)
for use in chemical knowledge graph construction. The preprocessing includes
renaming columns, consolidating multiple files, and saving the processed data
in a consistent format. This step ensures uniformity and ease of access for
subsequent data analysis and integration processes.

Classes:
    NodeDataProcessor: Handles preprocessing of assay, protein, gene, compound, experimental context, and assay endpoint data.

Example Usage:
    >>> processor = NodeDataProcessor(data_dir='path/to/data')
    >>> processor.preprocess_assays()
    >>> processor.preprocess_proteins()
    >>> processor.preprocess_genes()
    >>> processor.preprocess_compounds()
"""

import datetime
import glob
import os
import pandas as pd
import logging

import yaml

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
        preprocess_experimental_contexts(): Processes and stamps ExperimentalContext node table.
        preprocess_assay_endpoints(): Processes and stamps AssayEndpoint node table.
    """

    def __init__(self, data_dir: str, schema_path: str = "config/node_schema.yml"):
        """
        Initializes the NodeDataProcessor with a directory path to manage the data files.

        Args:
            data_dir (str): The directory where the node data files are stored.
        """
        self.data_dir = data_dir        
        self.schema_version = None
        self.schema = None
        if os.path.exists(schema_path):
            with open(schema_path, "r") as f:
                self.schema = yaml.safe_load(f)
                self.schema_version = self.schema.get("schema_version")

    def _stamp_provenance(self, df: pd.DataFrame) -> pd.DataFrame:
        df["ProcessedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["ProcessedBy"] = "ChemGraphBuilder.NodeDataProcessor/1.0"
        if self.schema_version:
            df["NodeSchemaVersion"] = self.schema_version
        return df


    def _validate_required_columns(self, df: pd.DataFrame, node_type: str, required_cols: list[str] | None = None) -> None:
        if required_cols is None and self.schema:
            # find node config by label or type
            for _, node_cfg in self.schema["nodes"].items():
                if node_cfg.get("label") == node_type or node_type in node_cfg.get("aliases", []):
                    required_cols = node_cfg.get("required_columns", [])
                    break
        required_cols = required_cols or []
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            logging.warning(
                "NodeDataProcessor: %s is missing required columns: %s",
                node_type,
                ", ".join(missing),
            )

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
        df = self._stamp_provenance(df)
        self._validate_required_columns(df, "Assay")
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
        df = self._stamp_provenance(df)
        self._validate_required_columns(
            df, "Protein")
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
        df = self._stamp_provenance(df)
        self._validate_required_columns(df, "Gene")
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
        df = self._stamp_provenance(df)
        self._validate_required_columns(df, "Compound")
        df.to_csv(f"{output_file.replace('.csv', '_Processed.csv')}", index=False)


    def preprocess_experimental_contexts(self):
        """
        Processes ExperimentalContext data derived from assay metadata and saves it
        as a processed node table under Data/Nodes.

        Expected input:
            Data/Nodes/ExperimentalContext_Properties.csv

        Output:
            Data/Nodes/ExperimentalContext_Properties_Processed.csv
        """
        in_path = f"{self.data_dir}/Nodes/ExperimentalContext_Properties.csv"
        if not os.path.exists(in_path):
            logging.warning("ExperimentalContext input not found at %s (skipping).", in_path)
            return None

        df = pd.read_csv(in_path, low_memory=False)

        # Type normalization
        for col in ("AssayID", "TaxonomyID"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

        df = self._stamp_provenance(df)
        self._validate_required_columns(df, "ExperimentalContext")

        out_path = f"{self.data_dir}/Nodes/ExperimentalContext_Properties_Processed.csv"
        df.to_csv(out_path, index=False)
        logging.info("Saved processed ExperimentalContext to %s (%d rows)", out_path, len(df))
        return df


    def preprocess_assay_endpoints(self):
        """
        Processes AssayEndpoint data derived from AllDataConnected assay results and
        saves it as a processed node table under Data/Nodes.

        Expected input:
            Data/Nodes/AssayEndpoint_Properties.csv

        Output:
            Data/Nodes/AssayEndpoint_Properties_Processed.csv
        """
        in_path = f"{self.data_dir}/Nodes/AssayEndpoint_Properties.csv"
        if not os.path.exists(in_path):
            logging.warning("AssayEndpoint input not found at %s (skipping).", in_path)
            return None

        df = pd.read_csv(in_path, low_memory=False)

        # Type normalization
        if "AssayID" in df.columns:
            df["AssayID"] = pd.to_numeric(df["AssayID"], errors="coerce").astype("Int64")

        if "HasNumericValue" in df.columns:
            # robust bool conversion
            df["HasNumericValue"] = df["HasNumericValue"].astype(str).str.lower().isin(["true", "1", "yes"])

        for col in ("ValueMin_uM", "ValueMedian_uM", "ValueMax_uM"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df = self._stamp_provenance(df)
        self._validate_required_columns(df, "AssayEndpoint")

        out_path = f"{self.data_dir}/Nodes/AssayEndpoint_Properties_Processed.csv"
        df.to_csv(out_path, index=False)
        logging.info("Saved processed AssayEndpoint to %s (%d rows)", out_path, len(df))
        return df
