"""
Module to collect and process relationship data for different types of relationships using RelationshipPropertiesExtractor and RelationshipDataProcessor.

Classes:
    RelationshipCollectorProcessor: Class to collect and process relationship data for different types of relationships.

Functions:
    main: Main function to parse command-line arguments and collect relationship data for the specified type.
"""

import os
import logging
import argparse
from relationship_properties_extractor import RelationshipPropertiesExtractor
from relationship_data_processor import RelationshipDataProcessor
from setup_data_folder import SetupDataFolder

logging.basicConfig(level=logging.INFO)

class RelationshipCollectorProcessor:
    """
    A class to collect and process relationship data for different types of relationships using RelationshipPropertiesExtractor and RelationshipDataProcessor.

    Attributes:
        relationship_type (str): The type of relationship to collect data for.
        data_file (str): The path to the data file containing relationship data.
        extractor (RelationshipPropertiesExtractor): An instance of RelationshipPropertiesExtractor.
        processor (RelationshipDataProcessor): An instance of RelationshipDataProcessor.
    """

    def __init__(self, relationship_type, start_chunk=0):
        """
        Initializes the RelationshipCollectorProcessor with the relationship type, data file, and start chunk index.

        Args:
            relationship_type (str): The type of relationship to collect data for (e.g., 'Assay_Compound', 'Assay_Enzyme', 'Gene_Enzyme', 'Compound_Enzyme', 'Compound_Similarity', 'Compound_Cooccurrence', 'Compound_Transformation').
            start_chunk (int): The starting chunk index for processing data.
        """
        self.relationship_type = relationship_type
        self.data_file = "Data/AllDataConnected.csv"
        self.start_chunk = start_chunk
        self.extractor = RelationshipPropertiesExtractor()
        self.processor = RelationshipDataProcessor(path="Data/Relationships/Assay_Compound_Relationship")

        # Setup data folder
        data_folder_setup = SetupDataFolder()
        data_folder_setup.setup()

    def collect_relationship_data(self):
        """
        Collects and processes relationship data based on the relationship type and saves it to the appropriate file.
        """
        if self.relationship_type == 'Assay_Compound':
            self.extractor.assay_compound_relationship(self.data_file, start_chunk=self.start_chunk)
            self.processor.process_files()
        elif self.relationship_type == 'Assay_Enzyme':
            self.extractor.assay_enzyme_relationship(self.data_file)
        elif self.relationship_type == 'Gene_Enzyme':
            self.extractor.gene_enzyme_relationship(self.data_file)
        elif self.relationship_type == 'Compound_Gene':
            self.extractor.compound_gene_relationship(self.data_file)
        elif self.relationship_type == 'Compound_Similarity':
            self.extractor.compound_similarity_relationship(self.data_file, start_chunk=self.start_chunk)
        elif self.relationship_type == 'Compound_Cooccurrence':
            self.extractor.compound_cooccurrence(self.data_file, start_chunk=self.start_chunk)
        elif self.relationship_type == 'Compound_Transformation':
            self.extractor.compound_transformation(self.data_file)
        else:
            logging.error(f"Unsupported relationship type: {self.relationship_type}")

def main():
    """
    Main function to parse command-line arguments and collect relationship data for the specified type.
    """
    parser = argparse.ArgumentParser(description="Collect relationship data for different types of relationships.")
    parser.add_argument('--relationship_type',
                        type=str,
                        required=True,
                        choices=['Assay_Compound', 'Assay_Enzyme', 'Gene_Enzyme',
                                 'Compound_Gene', 'Compound_Similarity',
                                 'Compound_Cooccurrence',
                                 'Compound_Transformation'],
                        help='The type of relationship to collect data for')
    parser.add_argument('--start_chunk',
                        type=int,
                        default=0,
                        help='The starting chunk index for processing data')

    args = parser.parse_args()

    collector = RelationshipCollectorProcessor(relationship_type=args.relationship_type, start_chunk=args.start_chunk)
    collector.collect_relationship_data()

if __name__ == '__main__':
    main()
