"""
Module to set up a data directory with a predefined structure.

This module provides the DataFolderSetup class, which creates a directory
structure for a data folder. The structure includes nodes and
relationships folders with specified subfolders.

Classes:
    DataFolderSetup: Class to set up a data directory with a predefined structure.

Functions:
    main: Main function to set up the data directory.
"""

import os

class SetupDataFolder:
    """
    Class to set up a data directory with a predefined structure.

    Attributes:
        data_folder (str): The name of the data folder.
        base_path (str): The base path for the data directory.
        structure (dict): The structure of directories to create.
    """

    def __init__(self):
        """
        Initializes the DataFolderSetup with the data folder name and directory structure.
        """
        self.data_folder = "Data"
        self.base_path = os.path.join(os.getcwd(), self.data_folder)
        self.structure = {
            "Nodes": ["Compound_properties"],
            "Relationships": [
                "Assay_Compound_Relationship",
                "Compound_Similarities",
                "Cpd_Cpd_CoOccurrence",
                "Cpd_Gene_CoOccurrence"
            ]
        }

    @staticmethod
    def create_folder(path):
        """
        Creates a folder if it does not already exist.

        Args:
            path (str): The path of the folder to create.
        """
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Created folder: {path}")
        else:
            print(f"Folder already exists: {path}")

    def setup(self):
        """
        Sets up the data directory structure based on the predefined structure.
        """
        # Create the base data directory
        self.create_folder(self.base_path)

        # Create the 'Nodes' directory and its subdirectories
        nodes_path = os.path.join(self.base_path, "Nodes")
        self.create_folder(nodes_path)
        for folder in self.structure["Nodes"]:
            self.create_folder(os.path.join(nodes_path, folder))

        # Create the 'Relationships' directory and its subdirectories
        relationships_path = os.path.join(self.base_path, "Relationships")
        self.create_folder(relationships_path)
        for folder in self.structure["Relationships"]:
            self.create_folder(os.path.join(relationships_path, folder))

        # Change the current working directory to the base data directory
        os.chdir(self.base_path)
        print(f"Changed current directory to: {self.base_path}")

def main():
    """
    Main function to set up the data directory.
    """
    data_folder_setup = SetupDataFolder()
    data_folder_setup.setup()

if __name__ == '__main__':
    main()
