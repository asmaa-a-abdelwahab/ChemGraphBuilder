"""
Module to load data into a Neo4j graph database for different node types.

This module provides the GraphDataLoader class, which allows loading data for
specific node types into a Neo4j database.
Users can provide the connection details and node label to load the data.

Classes:
    GraphDataLoader: Class to load data into a Neo4j graph database for different node types.

Functions:
    main: Main function to parse command-line arguments and load data for the specified node type.
"""

import logging
import argparse
import pandas as pd
import yaml
from neo4j import GraphDatabase
from chemgraphbuilder.add_graph_nodes import AddGraphNodes

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

class GraphNodesLoader:
    """
    Class to load data into a Neo4j graph database for different node types.

    Attributes:
        driver: The Neo4j driver instance.
        node_data_adder: An instance of the AddGraphNodes class.
        label_mapping: A dictionary mapping node labels to their unique properties and file paths.
    """

    def __init__(self, uri, username, password, schema_path="config/node_schema.yml"):
        """
        Initializes the GraphDataLoader with Neo4j connection details.

        Args:
            uri (str): The URI of the Neo4j database.
            username (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database.
        """
        self.driver = GraphDatabase.driver(uri, auth=(username, password))
        self.logger = logging.getLogger(__name__)  # Define the logger
        self.logger.info("GraphNodesLoader class initialized.")
        self.node_data_adder = AddGraphNodes(self.driver)
        with open(schema_path, "r") as f:
            schema = yaml.safe_load(f)
        # store only what we need
        self.label_mapping = {
            node_cfg["label"]: {
                "unique_property": node_cfg["unique_property"],
                "file_path": node_cfg["csv_ontologies"] or node_cfg["csv_processed"],
            }
            for _, node_cfg in schema["nodes"].items()
        }

    def validate_node_type(self, label: str) -> bool:
        """
        Validate that the CSV for a given node label:
          - exists
          - contains the unique_property
          - has no duplicate values in unique_property
        """
        if label not in self.label_mapping:
            self.logger.error("No mapping found for label: %s", label)
            return False

        unique_property = self.label_mapping[label]["unique_property"]
        file_path = self.label_mapping[label]["file_path"]

        try:
            df = pd.read_csv(file_path)
        except FileNotFoundError:
            self.logger.error("CSV file for label %s not found: %s", label, file_path)
            return False

        if unique_property not in df.columns:
            self.logger.error(
                "Unique property '%s' not found in %s for label %s",
                unique_property,
                file_path,
                label,
            )
            return False

        dup_count = df[unique_property].duplicated().sum()
        if dup_count > 0:
            self.logger.warning(
                "Found %d duplicate %s values in %s",
                dup_count,
                unique_property,
                file_path,
            )

        self.logger.info("Validation OK for label %s (file: %s)", label, file_path)
        return True

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

    def _upsert_node(self, tx, label, unique_property, props: dict):
        query = f"""
        MERGE (n:{label} {{ {unique_property}: $unique_id }})
        SET n += $props
        """
        unique_id = props.pop(unique_property)
        tx.run(query, unique_id=unique_id, props=props)
    
    def process_and_add_nodes(self, file_path, label, unique_property):
        df = pd.read_csv(file_path).dropna(subset=[unique_property], how="any")
        self.logger.info("Processing %d %s nodes", len(df), label)

        with self.driver.session() as session:
            def batch(tx):
                for _, row in df.iterrows():
                    self._upsert_node(
                        tx,
                        label=label,
                        unique_property=unique_property,
                        props=row.to_dict(),
                    )
            session.execute_write(batch)

    def load_data_for_node_type(self, label):
        """
        Loads data for a specific node type into the Neo4j database.

        Args:
            label (str): The label of the node.
        """
        if label not in self.label_mapping:
            self.logger.error("No mapping found for label: %s", label)
            return

        unique_property = self.label_mapping[label]["unique_property"]
        file_path = self.label_mapping[label]["file_path"]

        self.create_uniqueness_constraint(label, unique_property)
        self.process_and_add_nodes(file_path, label, unique_property)


def main():
    """
    Main function to parse command-line arguments and load data for the specified node type.
    """
    parser = argparse.ArgumentParser(description="Load data intoNeo4j graph database.")
    parser.add_argument("--uri", required=True, help="URI for the Neo4j database")
    parser.add_argument(
        "--username", required=True, help="Username for the Neo4j database"
    )
    parser.add_argument(
        "--password", required=True, help="Password for the Neo4j database"
    )
    parser.add_argument("--label", required=True, help="Label of the node")
    parser.add_argument(
        "--dry_run", action="store_true",
        help="If set, validate CSV and constraints but do not write nodes."
    )

    args = parser.parse_args()

    # Create an instance of GraphDataLoader and load data for the specified node type
    graph_nodes_loader = GraphNodesLoader(args.uri, args.username, args.password)
    if args.dry_run:
        graph_nodes_loader.validate_node_type(args.label)
    else:
        graph_nodes_loader.load_data_for_node_type(args.label)


if __name__ == "__main__":
    main()
