"""
Module to set up a project directory with a predefined structure.

This module provides the ProjectSetup class, which creates a directory
structure for a given project name. The structure includes nodes and
relationships folders with specified subfolders.

Classes:
    ProjectSetup: Class to set up a project directory with a predefined structure.

Functions:
    main: Main function to parse command-line arguments and set up the project directory.
"""

import os
import argparse

class ProjectSetup:
    """
    Class to set up a project directory with a predefined structure.
    
    Attributes:
        project_name (str): The name of the project.
        base_path (str): The base path for the project directory.
        structure (dict): The structure of directories to create.
    """

    def __init__(self, project_name):
        """
        Initializes the ProjectSetup with the project name and directory structure.
        
        Args:
            project_name (str): The name of the project.
        """
        self.project_name = project_name
        self.base_path = os.path.join(os.getcwd(), project_name)
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
        Sets up the project directory structure based on the predefined structure.
        """
        # Create the base project directory
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

        # Change the current working directory to the base project directory
        os.chdir(self.base_path)
        print(f"Changed current directory to: {self.base_path}")

def main():
    """
    Main function to parse command-line arguments and set up the project directory.
    """
    parser = argparse.ArgumentParser(description="Setup project directory structure.")
    parser.add_argument('project_name', type=str, help='Name of the project')
    args = parser.parse_args()

    project_setup = ProjectSetup(args.project_name)
    project_setup.setup()

if __name__ == '__main__':
    main()
