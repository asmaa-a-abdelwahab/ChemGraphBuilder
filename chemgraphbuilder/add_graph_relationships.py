"""
Module for adding relationship data from CSV files to a Neo4j database.

This module provides a class and methods to read relationship data from CSV files
and add them to a Neo4j database, including generating Cypher queries.
"""

import os
import glob
import ast
import json
import pandas as pd
import numpy as np
from chemgraphbuilder.neo4jdriver import Neo4jBase
import logging

# Set up logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class AddGraphRelationships(Neo4jBase):
    """
    A class used to add relationship data from a CSV file or a directory of CSV
    files to a Neo4j database.

    Methods
    -------
    generate_cypher_queries_from_file(file_path, rel_type, source_label, destination_label, rel_type_column=None):
        Generate Cypher queries to create relationships in Neo4j based on the data from the CSV file.
    generate_cypher_queries_from_directories(directory, rel_type, source_label, destination_label, rel_type_column=None):
        Generate Cypher queries to create relationships in Neo4j by merging CSV files from a directory.
    execute_queries(queries, batch_size=100):
        Execute a list of provided Cypher queries against the Neo4j database.
    combine_csv_files(input_directory):
        Combine multiple CSV files with the same columns into a single DataFrame.
    process_and_add_relationships(file_path, rel_type, source_label, destination_label, rel_type_column=None):
        Process the CSV file and add relationship data to the Neo4j database.
    process_and_add_relationships_from_directory(directory_path, rel_type, source_label, destination_label, rel_type_column=None):
        Combine CSV files from a directory and add relationship data to the Neo4j database.
    """

    def __init__(self, driver):
        """
        Initializes the AddGraphRelationships class with a Neo4j driver.

        Parameters
        ----------
        driver : neo4j.GraphDatabase.driver
            A driver instance to connect to the Neo4j database.
        """
        super().__init__()
        self.driver = driver
        self.logger.info("AddGraphRelationships class initialized.")


    @staticmethod
    def _generate_property_string(value):
        """
        Generate a property string for Cypher queries.

        Parameters
        ----------
        value : any
            The value to be converted to a string.

        Returns
        -------
        str
            The formatted string for the Cypher query.
        """
        if isinstance(value, (int, float)):
            return value
        try:
            return float(value)
        except (TypeError, ValueError):
            escaped_value = str(value).replace("\\", "\\\\").replace("'", "\\'").replace('"', '\\"').replace("\n", "\\\n")
            return f"'{escaped_value}'"


    def _process_properties(self, row, source_column, destination_column,
                            rel_type_column, standard_id):
        """
        Process the properties for the relationship from the CSV row.

        Parameters
        ----------
        row : pandas.Series
            The row from the DataFrame.
        source_column : str
            The column name for the source node.
        destination_column : str
            The column name for the destination node.
        rel_type_column : str
            The column name for the relationship type, if it is to be extracted from the CSV file.
        standard_id : dict
            A dictionary mapping standard IDs.

        Returns
        -------
        dict
            A dictionary of properties for the relationship.
        """
        columns_to_drop = [source_column, destination_column]
        if rel_type_column:
            columns_to_drop.append(rel_type_column)
        properties = row.drop(labels=columns_to_drop).to_dict()
        properties = {k: v for k, v in properties.items() if not pd.isna(v)}
        properties = {standard_id.get(k, k): v for k, v in properties.items()}
        properties = {f"`{k}`": v for k, v in properties.items()} ##last_edit
        return properties

    
    def _generate_query(self, source, target, relationship_type, properties,
                        source_label, destination_label, standard_id,
                        source_column, destination_column):
        """
        Generate a Cypher query for creating a relationship in Neo4j.

        Parameters
        ----------
        source : any
            The source node identifier.
        target : any
            The target node identifier.
        relationship_type : str
            The type of the relationship.
        properties : dict
            The properties of the relationship.
        source_label : str
            The label for the source node.
        destination_label : str
            The label for the destination node.
        standard_id : dict
            A dictionary mapping standard IDs.
        source_column : str
            The column name for the source node.
        destination_column : str
            The column name for the destination node.

        Returns
        -------
        str
            A Cypher query string.
        """
        source_value = self._generate_property_string(source)
        target_value = self._generate_property_string(target)
        source_key = standard_id[source_column] if source_column in standard_id else source_column
        destination_key = standard_id[destination_column] if destination_column in standard_id else destination_column

        query = (
            f"MATCH (a:{source_label} {{{source_key}: {source_value}}}), "
            f"(b:{destination_label} {{{destination_key}: {target_value}}}) "
            f"MERGE (a)-[r:{relationship_type}]->(b)"
        )

        if properties:
            set_clauses = [
                f"r.{prop.replace(' ', '').replace('[', '_').replace(']', '_')}"
                f"= {self._generate_property_string(value)}"
                for prop, value in properties.items()
            ]
            if set_clauses:
                query += " SET " + ", ".join(set_clauses)

        self.logger.debug(f"Generated query: {query}")
        return query


    def generate_cypher_queries_from_file(self, file_path, rel_type, source_label,
                                          destination_label, rel_type_column=None):
        """
        Generate Cypher queries for creating relationships in Neo4j based on the provided CSV file.

        Parameters
        ----------
        file_path : str
            The path to the CSV file.
        rel_type : str
            The type of the relationship.
        source_label : str
            The label for the source node.
        destination_label : str
            The label for the destination node.
        rel_type_column : str, optional
            The column name for the relationship type, if it is to be extracted from the CSV file.

        Yields
        ------
        str
            A Cypher query string.
        """
        standard_id = {
            'GeneID': 'GeneID',
            'Gene ID': 'GeneID',
            'Target GeneID': 'GeneID',
            'geneids': 'GeneID',
            'AssayID': 'AssayID',
            'Assay ID': 'AssayID',
            'AID': 'AssayID',
            'aid': 'AssayID',
            'ID_1': 'CompoundID',
            'ID_2': 'CompoundID',
            'substratecid': 'CompoundID',
            'metabolitecid': 'CompoundID',
            'Compound ID': 'CompoundID',
            'CompoundID': 'CompoundID',
            'CID': 'CompoundID',
            'Similar CIDs': 'CompoundID',
            'Target Accession': 'ProteinRefSeqAccession',
            'geneid': 'GeneID',
            'target_geneid': 'GeneID',
            'cid': 'CompoundID'
        }
        self.logger.info(f"Reading data from CSV file: {file_path}")
        # Read only the column names
        df = pd.read_csv(file_path, nrows=0)
        column_names = df.columns.tolist()
        # Step 1: Read the CSV file
        df = pd.read_csv(file_path, usecols=column_names, low_memory=False)
        # Step 2: Replace specific values with NaN
        df.replace("__nan__", np.nan, inplace=True)
        # Step 3: Drop columns that are completely empty (all NaN values)
        df = df.dropna(axis=1, how='all')
        df.to_csv("Data/test.csv", index=False)
        if rel_type == 'CO_OCCURS_IN_LITERATURE':
            df.rename(columns={df.columns[0]: list(ast.literal_eval(df[df.columns[0]][0]).keys())[0]},
                      inplace=True)
        source_column, destination_column = df.columns[:2] 
        df = df.dropna(subset=[source_column, destination_column], how='any')
        if df.empty:
            self.logger.error("The CSV file %s is empty or contains no valid data.",
                              file_path)
            return
        if rel_type_column:
            df = df.dropna(subset=[rel_type_column])

        for _, row in df.iterrows():
            source = row[source_column]
            destination = row[destination_column]
            relationship_type = row[rel_type_column] if rel_type_column else rel_type
            relationship_type = relationship_type.replace("/", "OR")
            properties = self._process_properties(row, source_column,
                                                  destination_column,
                                                  rel_type_column, standard_id)

            if rel_type == 'IS_SIMILAR_TO':
                targets = ast.literal_eval(row[destination_column])
                for target in targets:
                    query = self._generate_query(
                        source, target, relationship_type, properties,
                        source_label,destination_label, standard_id,
                        source_column, destination_column
                    )
                    yield query
                    
            elif rel_type == 'CO_OCCURS_IN_LITERATURE':
                source = ast.literal_eval(row[source_column])
                destination = ast.literal_eval(row[destination_column])
                if isinstance(source, dict):
                    source = list(source.values())[0]
                targets = ast.literal_eval(row[destination_column])
                if isinstance(targets, dict):
                    for target in targets.values():
                        query = self._generate_query(
                            source, target, relationship_type, properties,
                            source_label, destination_label, standard_id,
                            source_column, destination_column
                        )
                        yield query
                        
            elif rel_type == 'ENCODES':
                source = int(row[source_column])
                target = str(row[destination_column])
                query = self._generate_query(
                    source, target, relationship_type, properties, source_label,
                    destination_label, standard_id, source_column, destination_column
                )
                yield query
                
            else:
                source = int(row[source_column])
                target = int(row[destination_column])
                query = self._generate_query(
                    source, target, relationship_type, properties, source_label,
                    destination_label, standard_id, source_column, destination_column
                )
                yield query

        self.logger.info("Cypher queries generated successfully.")


    def generate_cypher_queries_from_directories(self, directory, rel_type,
                                                 source_label, destination_label,
                                                 rel_type_column=None):
        """
        Generate Cypher queries for creating relationships in Neo4j by processing each CSV file in a directory.
    
        Parameters
        ----------
        directory : str
            The path to the directory containing the CSV files.
        rel_type : str
            The type of the relationship.
        source_label : str
            The label for the source node.
        destination_label : str
            The label for the destination node.
        rel_type_column : str, optional
            The column name for the relationship type, if it is to be extracted from the CSV file.
    
        Yields
        ------
        str
            A Cypher query string.
        """    
        csv_files = glob.glob(directory)
        
        if not csv_files:
            self.logger.error("The directory %s contains no valid CSV files.", directory)
            return []
    
        for csv_file in csv_files:
            self.logger.info("Processing file: %s", csv_file)
            for query in self.generate_cypher_queries_from_file(csv_file, rel_type, source_label, destination_label, rel_type_column):
                yield query


    def execute_queries(self, queries, batch_size=100):
        """
        Execute the provided list of Cypher queries against the Neo4j database.

        Parameters
        ----------
        queries : list
            A list of Cypher query strings to execute.
        batch_size : int, optional
            The number of queries to execute in each batch.
        """
        if not queries:
            self.logger.error("No queries to execute.")
            return

        self.logger.info("Executing Cypher queries...")
        with self.driver.session() as session:
            for i in range(0, len(queries), batch_size):
                batch = queries[i:i + batch_size]
                try:
                    for query in batch:
                        session.run(query)
                        self.logger.debug(f"Executed query: {query}")
                    self.logger.info(f"Executed batch {i // batch_size + 1}"
                                     f"of {len(queries) // batch_size + 1}")
                except Exception as e:
                    self.logger.error("Failed to execute batch starting at query %s: %s",
                                      i, str(e))

        self.logger.info("All queries executed.")


    def process_and_add_relationships(self, file_path, rel_type, source_label,
                                      destination_label, rel_type_column=None):
        """
        Process the CSV file and add relationship data to the Neo4j database.

        Parameters
        ----------
        file_path : str
            The path to the CSV file.
        rel_type : str
            The type of the relationship.
        source_label : str
            The label for the source node.
        destination_label : str
            The label for the destination node.
        rel_type_column : str, optional
            The column name for the relationship type, if it is to be extracted from the CSV file.
        """
        self.logger.info("Processing and adding relationships from file: %s",
                         file_path)
        queries = list(self.generate_cypher_queries_from_file(file_path,
                                                              rel_type,
                                                              source_label,
                                                              destination_label,
                                                              rel_type_column))
        self.execute_queries(queries)
        self.logger.info("Successfully processed and added relationships from file: %s",
                         file_path)


    def process_and_add_relationships_from_directory(self, directory_path,
                                                     rel_type, source_label,
                                                     destination_label,
                                                     rel_type_column=None):
        """
        Combine CSV files from a directory and add relationship data to the Neo4j database.

        Parameters
        ----------
        directory_path : str
            The path to the directory containing the CSV files.
        rel_type : str
            The type of the relationship.
        source_label : str
            The label for the source node.
        destination_label : str
            The label for the destination node.
        rel_type_column : str, optional
            The column name for the relationship type, if it is to be extracted from the CSV file.
        """
        self.logger.info("Processing and adding relationships from directory: %s",
                         directory_path)
        queries = list(self.generate_cypher_queries_from_directories(directory_path,
                                                                     rel_type,
                                                                     source_label,
                                                                     destination_label,
                                                                     rel_type_column))
        if not queries:
            self.logger.error("No valid relationships found in the directory %s.",
                              directory_path)
            return

        self.execute_queries(queries)
        self.logger.info("Successfully processed and added relationships from directory: %s",
                         directory_path)

    def remove_self_relationships(self):
        """
        Remove self-relationships from the Neo4j database.
        """
        query = """
        MATCH (a)-[r]->(a)
        DELETE r
        """
        self.logger.info("Removing self-relationships from the database.")
        with self.driver.session() as session:
            try:
                session.run(query)
                self.logger.info("Self-relationships removed successfully.")
            except Exception as e:
                self.logger.error("Failed to remove self-relationships: %s", str(e))


    def make_relationships_bidirectional(self, relationship_type):
        """
        Ensure all relationships of a given type are bidirectional in the Neo4j database.

        Parameters
        ----------
        relationship_type : str
            The type of the relationship to make bidirectional.
        """
        query = f"""
        MATCH (a)-[r:{relationship_type}]->(b)
        WHERE NOT (b)-[:{relationship_type}]->(a)
        CREATE (b)-[:{relationship_type}]->(a)
        """
        self.logger.info(f"Making {relationship_type} relationships bidirectional.")
        with self.driver.session() as session:
            try:
                session.run(query)
                self.logger.info(f"Bidirectional relationships for {relationship_type} created successfully.")
            except Exception as e:
                self.logger.error(f"Failed to create bidirectional relationships for {relationship_type}: {str(e)}")
