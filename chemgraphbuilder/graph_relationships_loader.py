"""
GraphRelationshipsLoader class for loading graph relationships into a Neo4j database.
"""
import logging
import argparse
from neo4j import GraphDatabase
from chemgraphbuilder.add_graph_relationships import AddGraphRelationships

# Set up logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class GraphRelationshipsLoader:
    """
    Class for loading graph relationships into a Neo4j database using the AddGraphRelationships class.

    Attributes:
        uri (str): The URI for the Neo4j database.
        username (str): The username for the Neo4j database.
        password (str): The password for the Neo4j database.
        relationship_settings (dict): Predefined settings for different relationship types.
    """

    def __init__(self, uri, username, password):
        """
        Initializes the GraphRelationshipsLoader with Neo4j connection details.

        Args:
            uri (str): The URI for the Neo4j database.
            username (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database.
        """
        self.driver = GraphDatabase.driver(uri, auth=(username, password))
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.add_graph_relationships = AddGraphRelationships(self.driver)
        self.logger.info("GraphRelationshipsLoader class initialized.")

        # Predefined settings for different relationship types
        self.relationship_settings = {
            "Compound_Gene": {
                "file_path": "Data/Relationships/Compound_Gene_Relationship.csv",
                "source_label": "Compound",
                "destination_label": "Gene",
                "rel_type_column": "Activity"
            },
            "Assay_Compound": {
                "file_path": "Data/Relationships/Assay_Compound_Relationship.csv",
                "source_label": "BioAssay",
                "destination_label": "Compound",
                "rel_type_column": None,
                "relationship_type": "EVALUATES"
            },
            "Assay_Gene": {
                "file_path": "Data/Relationships/Assay_Gene_Relationship.csv",
                "source_label": "BioAssay",
                "destination_label": "Gene",
                "rel_type_column": None,
                "relationship_type": "STUDIES"
            },
            "Compound_Transformation": {
                "file_path": "Data/Relationships/Compound_Transformation.csv",
                "source_label": "Compound",
                "destination_label": "Compound",
                "rel_type_column": None,
                "relationship_type": "IS_METABOLIZED_TO"
            },
            "Gene_Enzyme": {
                "file_path": "Data/Relationships/Gene_Enzyme_Relationship.csv",
                "source_label": "Gene",
                "destination_label": "Protein",
                "rel_type_column": None,
                "relationship_type": "ENCODES"
            },
            "Compound_Similarities": {
                "file_path": "Data/Relationships/Compound_Similarities",
                "source_label": "Compound",
                "destination_label": "Compound",
                "rel_type_column": None,
                "relationship_type": "IS_SIMILAR_TO",
                "is_directory": True
            },
            "Cpd_Cpd_CoOccurence": {
                "file_path": "Data/Relationships/Cpd_Cpd_CoOccurence",
                "source_label": "Compound",
                "destination_label": "Compound",
                "rel_type_column": None,
                "relationship_type": "CO_OCCURS_IN_LITERATURE",
                "is_directory": True
            },
            "Cpd_Gene_CoOccurence": {
                "file_path": "Data/Relationships/Cpd_Gene_CoOccurence",
                "source_label": "Compound",
                "destination_label": "Gene",
                "rel_type_column": None,
                "relationship_type": "CO_OCCURS_IN_LITERATURE",
                "is_directory": True
            }
        }

    def close(self):
        """Closes the Neo4j database driver connection."""
        self.driver.close()

    def add_relationships(self, relationship_type):
        """
        Adds relationships to the Neo4j database based on the specified relationship type.

        Args:
            relationship_type (str): The type of the relationship.
        """
        settings = self.relationship_settings.get(relationship_type)
        if not settings:
            self.logger.error(f"Invalid relationship type: {relationship_type}")
            return

        file_path = settings["file_path"]
        source_label = settings["source_label"]
        destination_label = settings["destination_label"]
        rel_type_column = settings.get("rel_type_column")
        is_directory = settings.get("is_directory", False)

        if is_directory:
            self.add_graph_relationships.process_and_add_relationships_from_directory(file_path,
                                                                                      settings["relationship_type"],
                                                                                      source_label,
                                                                                      destination_label,
                                                                                      rel_type_column)
        else:
            self.add_graph_relationships.process_and_add_relationships(file_path,
                                                                       settings.get("relationship_type"),
                                                                       source_label,
                                                                       destination_label,
                                                                       rel_type_column)

def main():
    """Main function to handle command-line arguments and run the GraphRelationshipsLoader."""
    parser = argparse.ArgumentParser(description="Load graph relationships into"
                                                 "a Neo4j database.")
    parser.add_argument("--uri", type=str, required=True,
                        help="URI for the Neo4j database.")
    parser.add_argument("--username", type=str, required=True,
                        help="Username for the Neo4j database.")
    parser.add_argument("--password", type=str, required=True,
                        help="Password for the Neo4j database.")
    parser.add_argument("--relationship_type", type=str, required=True,
                        help="Type of the relationship to add.")
    args = parser.parse_args()

    loader = GraphRelationshipsLoader(args.uri, args.username, args.password)
    loader.add_relationships(args.relationship_type)
    loader.close()

if __name__ == "__main__":
    main()
