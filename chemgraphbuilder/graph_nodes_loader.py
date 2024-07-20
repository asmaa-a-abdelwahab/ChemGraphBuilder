"""
Module to load data into a Neo4j graph database for different node types.

This module provides the GraphDataLoader class, which allows loading data for specific node types into a Neo4j database.
Users can provide the connection details and node label to load the data.

Classes:
    GraphDataLoader: Class to load data into a Neo4j graph database for different node types.

Functions:
    main: Main function to parse command-line arguments and load data for the specified node type.
"""

import argparse
from neo4j import GraphDatabase
from chemgraphbuilder.add_graph_nodes import AddGraphNodes  # Assuming this is the module where AddGraphNodes is defined


class GraphNodesLoader:
    """
    Class to load data into a Neo4j graph database for different node types.

    Attributes:
        driver: The Neo4j driver instance.
        node_data_adder: An instance of the AddGraphNodes class.
        label_mapping: A dictionary mapping node labels to their unique properties and file paths.
    """

    def __init__(self, uri, username, password):
        """
        Initializes the GraphDataLoader with Neo4j connection details.

        Args:
            uri (str): The URI of the Neo4j database.
            username (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database.
        """
        self.driver = GraphDatabase.driver(uri, auth=(username, password))
        self.node_data_adder = AddGraphNodes(self.driver)
        self.label_mapping = {
            "Compound": {
                "unique_property": "CompoundID",
                "file_path": "Data/Nodes/Compound_Properties_Processed.csv"
            },
            "BioAssay": {
                "unique_property": "AssayID",
                "file_path": "Data/Nodes/Assay_Properties_Processed.csv"
            },
            "Gene": {
                "unique_property": "GeneID",
                "file_path": "Data/Nodes/Gene_Properties_Processed.csv"
            },
            "Protein": {
                "unique_property": "ProteinRefSeqAccession",
                "file_path": "Data/Nodes/Protein_Properties_Processed.csv"
            }
        }

    def create_uniqueness_constraint(self, label, unique_property):
        """
        Creates a uniqueness constraint for a given node label and property.

        Args:
            label (str): The label of the node.
            unique_property (str): The property to enforce uniqueness on.
        """
        self.node_data_adder.create_uniqueness_constraint(
            self.driver, label=label, unique_property=unique_property
        )

    def process_and_add_nodes(self, file_path, label, unique_property):
        """
        Processes and adds nodes from a CSV file to the Neo4j database.

        Args:
            file_path (str): The path to the CSV file containing node data.
            label (str): The label of the node.
            unique_property (str): The unique property of the node.
        """
        self.node_data_adder.process_and_add_nodes(
            file_path, label=label, unique_property=unique_property
        )

    def load_data_for_node_type(self, label):
        """
        Loads data for a specific node type into the Neo4j database.

        Args:
            label (str): The label of the node.
        """
        if label not in self.label_mapping:
            self.logger.error(f"No mapping found for label: {label}")
            return

        unique_property = self.label_mapping[label]["unique_property"]
        file_path = self.label_mapping[label]["file_path"]

        self.create_uniqueness_constraint(label, unique_property)
        self.process_and_add_nodes(file_path, label, unique_property)


def main():
    """
    Main function to parse command-line arguments and load data for the specified node type.
    """
    parser = argparse.ArgumentParser(description="Load data into Neo4j graph database.")
    parser.add_argument('--uri', required=True, help='URI for the Neo4j database')
    parser.add_argument('--username', required=True, help='Username for the Neo4j database')
    parser.add_argument('--password', required=True, help='Password for the Neo4j database')
    parser.add_argument('--label', required=True, help='Label of the node')

    args = parser.parse_args()

    # Create an instance of GraphDataLoader and load data for the specified node type
    graph_nodes_loader = GraphNodesLoader(args.uri, args.username, args.password)
    graph_nodes_loader.load_data_for_node_type(args.label)


if __name__ == "__main__":
    main()
