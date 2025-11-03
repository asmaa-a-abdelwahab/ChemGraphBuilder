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

import os
import logging
import argparse
from chemgraphbuilder.node_properties_extractor import NodePropertiesExtractor
from chemgraphbuilder.node_data_processor import NodeDataProcessor

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class NodesCollectorProcessor:
    """
    A class to collect and process data for different types of nodes using NodePropertiesExtractor and NodeDataProcessor.
    """

    def __init__(self, node_type, enzyme_list, start_chunk=None):
        """
        Initializes the NodesCollectorProcessor with the node type to collect data for, and the list of enzymes.

        Args:
            node_type (str): The type of node to collect data for (e.g., 'Compound', 'BioAssay', 'Gene', 'Protein').
            enzyme_list (list of str): List of enzyme names for which assay data will be fetched from PubChem.
            start_chunk (int, optional): The starting chunk index for processing. Default is None.
        """
        self.node_type = node_type
        self.extractor = NodePropertiesExtractor(enzyme_list=enzyme_list)
        self.processor = NodeDataProcessor(data_dir="Data")
        self.start_chunk = start_chunk

    def collect_and_process_data(self):
        """
        Collects and processes data based on the node type and saves it to the appropriate file.
        """
        data_file = "Data/AllDataConnected.csv"

        # Check if the data file exists before running the extractor
        if not os.path.exists(data_file):
            logging.info(f"{data_file} does not exist. Running data extraction...")
            df = self.extractor.run()
        else:
            logging.info(f"{data_file} already exists. Skipping main data extraction.")

        if self.node_type == "Compound":
            self.extractor.extract_compound_properties(
                main_data="Data/AllDataConnected.csv", start_chunk=self.start_chunk
            )
            self.processor.preprocess_compounds()
        elif self.node_type == "BioAssay":
            self.extractor.extract_assay_properties(
                main_data="Data/AllDataConnected.csv"
            )
            self.processor.preprocess_assays()
        elif self.node_type == "Gene":
            self.extractor.extract_gene_properties(
                main_data="Data/AllDataConnected.csv"
            )
            self.processor.preprocess_genes()
        elif self.node_type == "Protein":
            self.extractor.extract_protein_properties(
                main_data="Data/AllDataConnected.csv"
            )
            self.processor.preprocess_proteins()
        else:
            logging.error(f"Unsupported node type: {self.node_type}")


def main():
    """
    Main function to parse command-line arguments and collect data for the specified node type and enzyme list.
    """
    parser = argparse.ArgumentParser(
        description="Collect data for different types of nodes."
    )
    parser.add_argument(
        "--node_type",
        type=str,
        required=True,
        choices=["Compound", "BioAssay", "Gene", "Protein"],
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

    args = parser.parse_args()
    enzyme_list = args.enzyme_list.split(",")

    collector = NodesCollectorProcessor(
        node_type=args.node_type, enzyme_list=enzyme_list, start_chunk=args.start_chunk
    )
    collector.collect_and_process_data()


if __name__ == "__main__":
    main()