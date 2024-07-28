"""
NodeCollectorProcessor Module

This module provides the NodeCollectorProcessor class for collecting and processing data for different types of nodes
using the NodePropertiesExtractor and NodeDataProcessor classes. The collected data is intended for loading into
a Neo4j graph database. The module supports command-line interface (CLI) usage for ease of use.

Classes:
    NodeCollectorProcessor: A class to collect and process data for different types of nodes.

Functions:
    main: Main function to parse command-line arguments and collect data for the specified node type and enzyme list.
"""

import os
import logging
import argparse
from chemgraphbuilder.node_properties_extractor import NodePropertiesExtractor
from chemgraphbuilder.node_data_processor import NodeDataProcessor

logging.basicConfig(level=logging.INFO)

class NodeCollectorProcessor:
    """
    A class to collect and process data for different types of nodes using NodePropertiesExtractor and NodeDataProcessor.
    """

    def __init__(self, node_type, enzyme_list):
        """
        Initializes the NodesCollectorProcessor with Neo4j connection details, the node type to collect data for, and the list of enzymes.

        Args:
            uri (str): The URI for the Neo4j database.
            username (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database.
            node_type (str): The type of node to collect data for (e.g., 'Compound', 'BioAssay', 'Gene', 'Protein').
            enzyme_list (list of str): List of enzyme names for which assay data will be fetched from PubChem.
        """
        self.node_type = node_type
        self.extractor = NodePropertiesExtractor(enzyme_list=enzyme_list)
        self.processor = NodeDataProcessor(data_dir="Data")

    def collect_and_process_data(self):
        """
        Collects and processes data based on the node type and saves it to the appropriate file.
        """
        if self.node_type == 'Compound':
            df = self.extractor.run()
            self.extractor.extract_compound_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_compounds()
        elif self.node_type == 'BioAssay':
            df = self.extractor.run()
            self.extractor.extract_assay_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_assays()
        elif self.node_type == 'Gene':
            df = self.extractor.run()
            self.extractor.extract_gene_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_genes()
        elif self.node_type == 'Protein':
            df = self.extractor.run()
            self.extractor.extract_protein_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_proteins()
        else:
            logging.error(f"Unsupported node type: {self.node_type}")

def main():
    """
    Main function to parse command-line arguments and collect data for the specified node type and enzyme list.
    """
    parser = argparse.ArgumentParser(description="Collect data for different types of nodes.")
    parser.add_argument('--node_type', type=str, required=True, choices=['Compound', 'BioAssay', 'Gene', 'Protein'], help='The type of node to collect data for')
    parser.add_argument('--enzyme_list', type=str, required=True, help='Comma-separated list of enzyme names to fetch data for')

    args = parser.parse_args()
    enzyme_list = args.enzyme_list.split(',')

    collector = NodeCollectorProcessor(node_type=args.node_type, enzyme_list=enzyme_list)
    collector.collect_and_process_data()

if __name__ == '__main__':
    main()
