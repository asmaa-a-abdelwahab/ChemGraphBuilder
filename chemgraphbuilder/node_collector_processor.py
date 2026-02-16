"""
NodesCollectorProcessor Module

This module provides the NodesCollectorProcessor class for collecting and processing data for different types of nodes
using the NodePropertiesExtractor and NodeDataProcessor classes. The collected data is intended for loading into
a Neo4j graph database. The module supports command-line interface (CLI) usage for ease of use.

Classes:
    NodesCollectorProcessor: A class to collect and process data for different types of nodes.

Functions:
    main: Main function to parse command-line arguments and collect data for the specified node type and enzyme list.
"""

import json
import os
import logging
import argparse

from chemgraphbuilder.node_properties_extractor import NodePropertiesExtractor
from chemgraphbuilder.node_data_processor import NodeDataProcessor
from chemgraphbuilder.node_ontology_enricher import NodesOntologyEnricher

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)



class NodesCollectorProcessor:
    """
    A class to collect and process data for different types of nodes using
    NodePropertiesExtractor, NodeDataProcessor, and (optionally) NodesOntologyEnricher.
    """

    def __init__(
        self,
        node_type,
        enzyme_list,
        start_chunk=None,
        with_ontologies: bool = False,
        species: str = "human",
        mygene_email: str | None = None,
    ):
        """
        Initializes the NodesCollectorProcessor.

        Args:
            node_type (str): 'Compound', 'BioAssay', 'Gene', 'Protein', 'ExperimentalContext', or 'AssayEndpoint'.
            enzyme_list (list[str]): Enzyme names to fetch data for.
            start_chunk (int, optional): Starting chunk index for compounds.
            with_ontologies (bool): If True, run ontology enrichment after processing.
            species (str): Species passed to ontology services (e.g. MyGene.info).
            mygene_email (str, optional): Contact email for MyGene.info client.
        """
        self.node_type = node_type
        self.extractor = NodePropertiesExtractor(enzyme_list=enzyme_list)
        self.processor = NodeDataProcessor(data_dir="Data")
        self.start_chunk = start_chunk
        self.with_ontologies = with_ontologies

        # Enriched tables are written under Data/Nodes by NodeDataProcessor
        self.ontology_enricher = None
        if self.with_ontologies:
            self.ontology_enricher = NodesOntologyEnricher(
                data_dir="Data/Nodes",
                species=species,
                mygene_email=mygene_email,
            )

    def collect_and_process_data(self) -> None:
        """
        Collect and process data for the configured node type, and optionally enrich it
        with ontology annotations.

        Steps:
            1. Ensure the main joined file (Data/AllDataConnected.csv) exists.
            2. For the selected node_type, run:
               - extraction (NodePropertiesExtractor)
               - preprocessing (NodeDataProcessor)
               - optional ontology enrichment (NodesOntologyEnricher)
            3. Log a lightweight summary (row counts and file paths).
        """
        data_file = "Data/AllDataConnected.csv"

        # Helper: safe CSV row counter
        def _safe_count_csv(path: str) -> int:
            if not os.path.exists(path):
                return 0
            try:
                import pandas as pd  # local import to avoid hard dependency if unused
                return len(pd.read_csv(path))
            except Exception as exc:  # pragma: no cover - defensive
                logging.warning("Could not count rows in %s: %s", path, exc)
                return 0

        # 1. Ensure the main joined file exists
        if not os.path.exists(data_file):
            logging.info("%s does not exist. Running main data extraction...", data_file)
            # This will fetch assays + linked genes/proteins/compounds for the enzyme_list
            self.extractor.run()
        else:
            logging.info("%s already exists. Skipping main data extraction.", data_file)

        summary: dict[str, object] = {
            "node_type": self.node_type,
            "main_data_file": data_file,
        }

        # ------------------------------------------------------------------
        # 2. Branch per node type
        # ------------------------------------------------------------------

        # ---------------------- COMPOUNDS ----------------------
        if self.node_type == "Compound":
            logging.info("Starting COMPOUND extraction & preprocessing...")
            self.extractor.extract_compound_properties(
                main_data=data_file,
                start_chunk=self.start_chunk,
            )

            # Preprocess compound node table(s)
            self.processor.preprocess_compounds()
            processed_path = "Data/Nodes/Compound_Properties_Processed.csv"
            summary["processed_file"] = processed_path
            summary["n_processed_rows"] = _safe_count_csv(processed_path)

            # Ontology enrichment for compound nodes
            if self.ontology_enricher is not None:
                logging.info("Running compound ontology enrichment (ChEBI/MeSH/ChEMBL)...")
                # IMPORTANT: use the processed column name 'CompoundID'
                df_enriched = self.ontology_enricher.enrich_compounds(
                    input_name=os.path.basename(processed_path),
                    output_name="Compound_Properties_WithOntologies.csv",
                )
                if df_enriched is not None:
                    enriched_path = "Data/Nodes/Compound_Properties_WithOntologies.csv"
                    summary["enriched_file"] = enriched_path
                    summary["n_enriched_rows"] = len(df_enriched)

        # ---------------------- BIOASSAYS ----------------------
        elif self.node_type == "BioAssay":
            logging.info("Starting BIOASSAY extraction & preprocessing...")
            self.extractor.extract_assay_properties(main_data=data_file)

            self.processor.preprocess_assays()
            processed_path = "Data/Nodes/Assay_Properties_Processed.csv"
            summary["processed_file"] = processed_path
            summary["n_processed_rows"] = _safe_count_csv(processed_path)

            # Ontology enrichment for assay nodes (BAO)
            if self.ontology_enricher is not None:
                logging.info("Running assay ontology enrichment (BAO via OLS)...")
                # IMPORTANT: use the *processed* column names from NodeDataProcessor
                df_enriched = self.ontology_enricher.enrich_assays(
                    input_name=os.path.basename(processed_path),
                    output_name="Assay_Properties_WithOntologies.csv",
                    name_cols=[
                        "AssayName",
                        "AssayActivityName",
                        "AssayType",
                        "AssayDescription",
                    ],
                )
                if df_enriched is not None:
                    enriched_path = "Data/Nodes/Assay_Properties_WithOntologies.csv"
                    summary["enriched_file"] = enriched_path
                    summary["n_enriched_rows"] = len(df_enriched)

            # ------------------------------------------------------------------
            # Derived schema-support nodes (ExperimentalContext + AssayEndpoint)
            # These are computed locally from AllDataConnected + processed assay table.
            # ------------------------------------------------------------------
            try:
                logging.info("Deriving ExperimentalContext and AssayEndpoint node tables...")
                self.extractor.extract_experimental_context_properties(
                    main_data=data_file,
                    assay_properties_processed="Data/Nodes/Assay_Properties_Processed.csv",
                    gene_properties_processed="Data/Nodes/Gene_Properties_Processed.csv",
                    out_path="Data/Nodes/ExperimentalContext_Properties.csv",
                )
                self.processor.preprocess_experimental_contexts()

                self.extractor.extract_assay_endpoint_properties(
                    main_data=data_file,
                    assay_properties_processed="Data/Nodes/Assay_Properties_Processed.csv",
                    out_path="Data/Nodes/AssayEndpoint_Properties.csv",
                )
                self.processor.preprocess_assay_endpoints()


                # Record derived node outputs in the parent BioAssay summary
                ec_processed = "Data/Nodes/ExperimentalContext_Properties_Processed.csv"
                aep_processed = "Data/Nodes/AssayEndpoint_Properties_Processed.csv"
                summary["derived_nodes"] = {
                    "ExperimentalContext": {
                        "node_type": "ExperimentalContext",
                        "processed_file": ec_processed,
                        "n_processed_rows": _safe_count_csv(ec_processed),
                    },
                    "AssayEndpoint": {
                        "node_type": "AssayEndpoint",
                        "processed_file": aep_processed,
                        "n_processed_rows": _safe_count_csv(aep_processed),
                    },
                }

                # Optional ontology enrichment if supported by NodesOntologyEnricher implementation
                if self.ontology_enricher is not None:
                    if hasattr(self.ontology_enricher, "enrich_experimental_contexts"):
                        self.ontology_enricher.enrich_experimental_contexts(
                            input_name="ExperimentalContext_Properties_Processed.csv",
                            output_name="ExperimentalContext_Properties_WithOntologies.csv",
                        )
                    if hasattr(self.ontology_enricher, "enrich_assay_endpoints"):
                        self.ontology_enricher.enrich_assay_endpoints(
                            input_name="AssayEndpoint_Properties_Processed.csv",
                            output_name="AssayEndpoint_Properties_WithOntologies.csv",
                        )
            except Exception as exc:  # pragma: no cover
                logging.warning("Could not derive ExperimentalContext/AssayEndpoint: %s", exc)

        # ---------------------- GENES ----------------------
        elif self.node_type == "Gene":
            logging.info("Starting GENE extraction & preprocessing...")
            self.extractor.extract_gene_properties(main_data=data_file)

            self.processor.preprocess_genes()
            processed_path = "Data/Nodes/Gene_Properties_Processed.csv"
            summary["processed_file"] = processed_path
            summary["n_processed_rows"] = _safe_count_csv(processed_path)

            # Ontology enrichment for gene nodes (HGNC, Ensembl, GO, UniProt)
            if self.ontology_enricher is not None:
                logging.info("Running gene ontology enrichment (MyGene.info)...")
                df_enriched = self.ontology_enricher.enrich_genes(
                    input_name=os.path.basename(processed_path),
                    output_name="Gene_Properties_WithOntologies.csv",
                    gene_id_col="GeneID",
                )
                if df_enriched is not None:
                    enriched_path = "Data/Nodes/Gene_Properties_WithOntologies.csv"
                    summary["enriched_file"] = enriched_path
                    summary["n_enriched_rows"] = len(df_enriched)

        # ---------------------- PROTEINS ----------------------
        elif self.node_type == "Protein":
            logging.info("Starting PROTEIN extraction & preprocessing...")
            self.extractor.extract_protein_properties(main_data=data_file)

            self.processor.preprocess_proteins()
            processed_path = "Data/Nodes/Protein_Properties_Processed.csv"
            summary["processed_file"] = processed_path
            summary["n_processed_rows"] = _safe_count_csv(processed_path)

            # Ontology enrichment for protein nodes (UniProt, GO, InterPro)
            if self.ontology_enricher is not None:
                logging.info("Running protein ontology enrichment (UniProt/GO/InterPro)...")
                # IMPORTANT: use the processed column name 'ProteinRefSeqAccession'
                df_enriched = self.ontology_enricher.enrich_proteins(
                    input_name=os.path.basename(processed_path),
                    output_name="Protein_Properties_WithOntologies.csv",
                    accession_col="ProteinRefSeqAccession",
                )
                if df_enriched is not None:
                    enriched_path = "Data/Nodes/Protein_Properties_WithOntologies.csv"
                    summary["enriched_file"] = enriched_path
                    summary["n_enriched_rows"] = len(df_enriched)

        # ---------------------- EXPERIMENTAL CONTEXT ----------------------
        elif self.node_type == "ExperimentalContext":
            logging.info("Starting EXPERIMENTAL CONTEXT derivation & preprocessing...")
            self.extractor.extract_experimental_context_properties(
                main_data=data_file,
                assay_properties_processed="Data/Nodes/Assay_Properties_Processed.csv",
                gene_properties_processed="Data/Nodes/Gene_Properties_Processed.csv",
                out_path="Data/Nodes/ExperimentalContext_Properties.csv",
            )
            self.processor.preprocess_experimental_contexts()
            processed_path = "Data/Nodes/ExperimentalContext_Properties_Processed.csv"
            summary["processed_file"] = processed_path
            summary["n_processed_rows"] = _safe_count_csv(processed_path)

            if self.ontology_enricher is not None and hasattr(self.ontology_enricher, "enrich_experimental_contexts"):
                df_enriched = self.ontology_enricher.enrich_experimental_contexts(
                    input_name=os.path.basename(processed_path),
                    output_name="ExperimentalContext_Properties_WithOntologies.csv",
                )
                if df_enriched is not None:
                    enriched_path = "Data/Nodes/ExperimentalContext_Properties_WithOntologies.csv"
                    summary["enriched_file"] = enriched_path
                    summary["n_enriched_rows"] = len(df_enriched)

        # ---------------------- ASSAY ENDPOINT ----------------------
        elif self.node_type == "AssayEndpoint":
            logging.info("Starting ASSAY ENDPOINT derivation & preprocessing...")
            self.extractor.extract_assay_endpoint_properties(
                main_data=data_file,
                assay_properties_processed="Data/Nodes/Assay_Properties_Processed.csv",
                out_path="Data/Nodes/AssayEndpoint_Properties.csv",
            )
            self.processor.preprocess_assay_endpoints()
            processed_path = "Data/Nodes/AssayEndpoint_Properties_Processed.csv"
            summary["processed_file"] = processed_path
            summary["n_processed_rows"] = _safe_count_csv(processed_path)

            if self.ontology_enricher is not None and hasattr(self.ontology_enricher, "enrich_assay_endpoints"):
                df_enriched = self.ontology_enricher.enrich_assay_endpoints(
                    input_name=os.path.basename(processed_path),
                    output_name="AssayEndpoint_Properties_WithOntologies.csv",
                )
                if df_enriched is not None:
                    enriched_path = "Data/Nodes/AssayEndpoint_Properties_WithOntologies.csv"
                    summary["enriched_file"] = enriched_path
                    summary["n_enriched_rows"] = len(df_enriched)

        else:
            raise ValueError(
                f"Unsupported node_type '{self.node_type}'. "
                "Expected one of: 'Compound', 'BioAssay', 'Gene', 'Protein', 'ExperimentalContext', 'AssayEndpoint'."
            )

        # ------------------------------------------------------------------
        # Persist run summary JSON(s)
        # ------------------------------------------------------------------
        report_path = f"Data/Nodes/{self.node_type}_run_summary.json"
        try:
            with open(report_path, "w", encoding="utf-8") as f:
                json.dump(summary, f, indent=2)
            logging.info("Saved node run summary to %s", report_path)
        except Exception as exc:  # pragma: no cover
            logging.warning("Could not write run summary to %s: %s", report_path, exc)

        # If BioAssay derived schema-support nodes, write their own summaries too (if present)
        derived = summary.get("derived_nodes")
        if isinstance(derived, dict):
            for derived_type, derived_summary in derived.items():
                if not isinstance(derived_summary, dict):
                    continue
                dpath = f"Data/Nodes/{derived_type}_run_summary.json"
                try:
                    with open(dpath, "w", encoding="utf-8") as f:
                        json.dump(derived_summary, f, indent=2)
                    logging.info("Saved derived node run summary to %s", dpath)
                except Exception as exc:  # pragma: no cover
                    logging.warning("Could not write derived run summary to %s: %s", dpath, exc)



def main():
    """
    Main function to parse command-line arguments and collect data for the
    specified node type and enzyme list.
    """
    parser = argparse.ArgumentParser(
        description="Collect data for different types of nodes."
    )
    parser.add_argument(
        "--node_type",
        type=str,
        required=True,
        choices=["Compound", "BioAssay", "Gene", "Protein", "ExperimentalContext", "AssayEndpoint"],
        help="The type of node to collect data for",
    )
    parser.add_argument(
        "--enzyme_list",
        type=str,
        required=True,
        help="Comma-separated list of enzyme names to fetch data for",
    )
    parser.add_argument(
        "--start_chunk",
        type=int,
        default=None,
        help="The starting chunk index for processing Compound Data",
    )
    parser.add_argument(
        "--with_ontologies",
        action="store_true",
        help="If set, enrich processed node tables with ontology IDs "
             "(ChEBI / HGNC / Ensembl / UniProt / GO / BAO, etc.)",
    )
    parser.add_argument(
        "--species",
        type=str,
        default="human",
        help="Species passed to ontology services (e.g. 'human', 'mouse').",
    )
    parser.add_argument(
        "--mygene_email",
        type=str,
        default=None,
        help="Optional email used for MyGene.info client (polite, not required).",
    )

    args = parser.parse_args()
    enzyme_list = args.enzyme_list.split(",")

    collector = NodesCollectorProcessor(
        node_type=args.node_type,
        enzyme_list=enzyme_list,
        start_chunk=args.start_chunk,
        with_ontologies=args.with_ontologies,
        species=args.species,
        mygene_email=args.mygene_email,
    )
    collector.collect_and_process_data()

if __name__ == "__main__":
    main()