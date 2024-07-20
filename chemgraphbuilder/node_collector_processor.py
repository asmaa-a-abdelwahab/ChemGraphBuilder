"""
NodesCollectorProcessor Module

This module provides the NodesCollectorProcessor class for collecting and processing data for different types of nodes 
using the NodePropertiesExtractor and NodeDataProcessor classes. The collected data is intended for loading into 
a Neo4j graph database. The module supports command-line interface (CLI) usage for ease of use.

Classes:
    NodesCollectorProcessor: A class to collect and process data for different types of nodes.
    
Functions:
    main: Main function to parse command-line arguments and collect data for the specified node type and enzyme list.

Usage Example:
    To use the module from the command line:
    ```
    python graph_data_collector.py --uri bolt://localhost:7689 --username neo4j --password your_password --node_type Compound --enzyme_list CYP2D6,CYP3A4
    ```

    To use the class in a Python script or Jupyter notebook:
    ```
    enzyme_list = ['CYP2D6', 'CYP3A4']
    collector = NodesCollectorProcessor(uri="bolt://localhost:7689", username="neo4j", password="your_password", node_type="Compound", enzyme_list=enzyme_list)
    collector.collect_and_process_data()
    collector.close()
    ```
"""

import os
import logging
import argparse
from neo4j import GraphDatabase
from chemgraphbuilder.node_properties_extractor import NodePropertiesExtractor
from chemgraphbuilder.node_data_processor import NodeDataProcessor

logging.basicConfig(level=logging.INFO)

class NodeCollectorProcessor:
    """
    A class to collect and process data for different types of nodes using NodePropertiesExtractor and NodeDataProcessor.
    """

    def __init__(self, uri, username, password, node_type, enzyme_list):
        """
        Initializes the NodesCollectorProcessor with Neo4j connection details, the node type to collect data for, and the list of enzymes.

        Args:
            uri (str): The URI for the Neo4j database.
            username (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database.
            node_type (str): The type of node to collect data for (e.g., 'Compound', 'BioAssay', 'Gene', 'Protein').
            enzyme_list (list of str): List of enzyme names for which assay data will be fetched from PubChem.
        """
        self.driver = GraphDatabase.driver(uri, auth=(username, password))
        self.node_type = node_type
        self.extractor = NodePropertiesExtractor(enzyme_list=enzyme_list)
        self.processor = NodeDataProcessor(data_dir="Data")

    def close(self):
        """Closes the connection to the Neo4j database."""
        self.driver.close()

    def collect_and_process_data(self):
        """
        Collects and processes data based on the node type and saves it to the appropriate file.
        """
        if self.node_type == 'Compound':
            self.extractor.create_data_directories()
            df = self.extractor.run()
            self.extractor.extract_compound_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_compounds()
        elif self.node_type == 'BioAssay':
            self.extractor.create_data_directories()
            df = self.extractor.run()
            self.extractor.extract_assay_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_assays()
        elif self.node_type == 'Gene':
            self.extractor.create_data_directories()
            df = self.extractor.run()
            self.extractor.extract_gene_properties(main_data='Data/AllDataConnected.csv')
            self.processor.preprocess_genes()
        elif self.node_type == 'Protein':
            self.extractor.create_data_directories()
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
    parser.add_argument('--uri', type=str, required=True, help='The URI for the Neo4j database')
    parser.add_argument('--username', type=str, required=True, help='The username for the Neo4j database')
    parser.add_argument('--password', type=str, required=True, help='The password for the Neo4j database')
    parser.add_argument('--node_type', type=str, required=True, choices=['Compound', 'BioAssay', 'Gene', 'Protein'], help='The type of node to collect data for')
    parser.add_argument('--enzyme_list', type=str, required=True, help='Comma-separated list of enzyme names to fetch data for')

    args = parser.parse_args()
    enzyme_list = args.enzyme_list.split(',')

    collector = NodesCollectorProcessor(uri=args.uri, username=args.username, password=args.password, node_type=args.node_type, enzyme_list=enzyme_list)
    collector.collect_and_process_data()
    collector.close()

if __name__ == '__main__':
    main()
