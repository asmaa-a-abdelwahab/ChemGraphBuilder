"""
Module for adding node data from CSV files to a Neo4j database.

This module provides a class and methods to read node data from CSV files
and add them to a Neo4j database, including creating uniqueness constraints
and generating Cypher queries.
"""

import os
import logging
import pandas as pd
from neo4jdriver import Neo4jBase

class AddGraphNodes(Neo4jBase):
    """
    A class used to add node data from a CSV file or a directory of CSV files to a Neo4j database.

    Methods:
    --------
    create_uniqueness_constraint(driver, label, unique_property):
        Create a uniqueness constraint for the unique property of nodes in Neo4j.
    generate_cypher_queries(node_dict, label, unique_property):
        Generate Cypher queries to update nodes in Neo4j based on the data from the CSV file.
    execute_queries(queries):
        Execute a list of provided Cypher queries against the Neo4j database.
    read_csv_file(file_path, unique_property):
        Read data from a CSV file and extract node properties.
    combine_csv_files(input_directory):
        Combine multiple CSV files with the same columns into a single DataFrame.
    process_and_add_nodes(file_path, label, unique_property):
        Process the CSV file and add node data to the Neo4j database.
    process_and_add_nodes_from_directory(directory_path, label, unique_property):
        Combine CSV files from a directory and add node data to the Neo4j database.
    """

    def __init__(self, driver):
        """
        Initializes the AddGraphNodes class with a Neo4j driver.

        Parameters:
        -----------
        driver : neo4j.GraphDatabase.driver
            A driver instance to connect to the Neo4j database.
        """
        super().__init__()
        self.driver = driver
        self.logger.info("AddGraphNodes class initialized.")

    @staticmethod
    def create_uniqueness_constraint(driver, label, unique_property):
        """
        Create a uniqueness constraint for the unique property of nodes in Neo4j.

        Parameters:
        -----------
        driver : neo4j.GraphDatabase.driver
            A driver instance to connect to the Neo4j database.
        label : str
            The label of the node.
        unique_property : str
            The unique property of the node.
        """
        constraint_query = (
            f"CREATE CONSTRAINT IF NOT EXISTS FOR (n:{label}) "
            f"REQUIRE n.{unique_property} IS UNIQUE"
        )
        with driver.session() as session:
            try:
                session.run(constraint_query)
                logging.info(
                    "Uniqueness constraint created successfully on %s property of %s nodes.",
                    unique_property, label)
            except Exception as e:
                logging.error("Failed to create uniqueness constraint: %s", e)

    @staticmethod
    def _generate_property_string(value):
        if isinstance(value, (int, float)):
            return value
        try:
            return float(value)
        except (TypeError, ValueError):
            escaped_value = "'" + str(value).replace("'", "\\'").replace("\n", "\\n") + "'"
            return escaped_value

    def generate_cypher_queries(self, node_dict, label, unique_property):
        """
        Generate Cypher queries for updating Neo4j based on the provided node data dictionary.

        Parameters:
        -----------
        node_dict : dict
            A dictionary with unique identifiers as keys and node data as values.
        label : str
            The label of the node.
        unique_property : str
            The unique property of the node.

        Yields:
        -------
        str
            A Cypher query string.
        """
        # Create an index for the unique_property
        create_index_query = f"CREATE INDEX IF NOT EXISTS FOR (n:{label}) ON (n.{unique_property})"
        self.logger.debug(create_index_query)
        yield create_index_query

        for unique_id, properties in node_dict.items():
            unique_id = f'"{unique_id}"' if isinstance(unique_id, str) else unique_id
            query = f"MERGE (n:{label} {{{unique_property}: {unique_id}}})"
            set_clauses = [
                f"n.{prop.replace(' ', '')} = {self._generate_property_string(value)}"
                for prop, value in properties.items()
            ]
            if set_clauses:
                query += " SET " + ", ".join(set_clauses)
            else:
                query += ";"
            self.logger.debug(query)
            yield query
        self.logger.info("Cypher queries generated successfully.")

    def execute_queries(self, queries):
        """
        Execute the provided list of Cypher queries against the Neo4j database.

        Parameters:
        -----------
        queries : list
            A list of Cypher query strings to execute.
        """
        self.logger.info("Executing Cypher queries...")
        with self.driver.session() as session:
            self.logger.info("Executing Cypher queries Started....")
            for query in queries:
                try:
                    session.run(query)
                except Exception as e:
                    self.logger.error("Failed to execute query: %s", e)
        self.logger.info("All queries executed.")

    def read_csv_file(self, file_path, unique_property):
        """
        Read data from a CSV file and extract node properties.

        Parameters:
        -----------
        file_path : str
            The path to the CSV file.
        unique_property : str
            The column name that serves as the unique identifier for the nodes.

        Returns:
        --------
        dict
            A dictionary with unique identifiers as keys and extracted data as values.
        """
        self.logger.info("Reading data from CSV file: %s", file_path)
        df = pd.read_csv(file_path).dropna(subset=[unique_property], how='any')
        node_dict = {
            row[unique_property]: row.drop(labels=[unique_property]).to_dict()
            for _, row in df.iterrows()
        }
        self.logger.info("Successfully read data for %d nodes from CSV.", len(node_dict))
        return node_dict

    def combine_csv_files(self, input_directory):
        """
        Combine multiple CSV files with the same columns into a single DataFrame.

        Parameters:
        -----------
        input_directory : str
            The directory containing the CSV files to be combined.

        Returns:
        --------
        DataFrame
            A combined DataFrame containing data from all the CSV files.
        """
        self.logger.info("Combining CSV files from directory: %s", input_directory)
        dfs = [
            pd.read_csv(os.path.join(input_directory, file))
            for file in os.listdir(input_directory)
            if file.endswith(".csv")
        ]
        combined_df = pd.concat(dfs, ignore_index=True)
        self.logger.info("Successfully combined %d CSV files.", len(dfs))
        return combined_df

    def process_and_add_nodes(self, file_path, label, unique_property):
        """
        Process the CSV file and add node data to the Neo4j database.

        Parameters:
        -----------
        file_path : str
            The path to the CSV file.
        label : str
            The label of the node.
        unique_property : str
            The unique property of the node.
        """
        self.logger.info("Processing and adding nodes from file: %s", file_path)
        node_dict = self.read_csv_file(file_path, unique_property)
        queries = list(self.generate_cypher_queries(node_dict, label, unique_property))
        self.execute_queries(queries)
        self.logger.info("Successfully processed and added nodes from file: %s", file_path)

    def process_and_add_nodes_from_directory(self, directory_path, label, unique_property):
        """
        Combine CSV files from a directory and add node data to the Neo4j database.

        Parameters:
        -----------
        directory_path : str
            The path to the directory containing the CSV files.
        label : str
            The label of the node.
        unique_property : str
            The unique property of the node.
        """
        self.logger.info("Processing and adding nodes from directory: %s", directory_path)
        combined_df = self.combine_csv_files(directory_path)
        temp_file = os.path.join(directory_path, "combined_temp.csv")
        combined_df.to_csv(temp_file, index=False)
        self.process_and_add_nodes(temp_file, label, unique_property)
        os.remove(temp_file)
        self.logger.info("Successfully processed and added nodes from directory: %s",
                         directory_path)

    def public_generate_property_string(self, value):
        """
        Public method to access the protected _generate_property_string method for testing.

        Parameters:
        -----------
        value : Any
            The value to be formatted.

        Returns:
        --------
        str
            The formatted property string.
        """
        return self._generate_property_string(value)
