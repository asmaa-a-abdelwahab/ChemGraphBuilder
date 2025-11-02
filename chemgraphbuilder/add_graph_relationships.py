"""
Module for adding relationship data from CSV files to a Neo4j database.

This module provides a class and methods to read relationship data from CSV files
and add them to a Neo4j database, including generating Cypher queries.
"""

import glob
import ast
import pandas as pd
import numpy as np
from chemgraphbuilder.neo4jdriver import Neo4jBase
import logging

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


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
        
    def _normalize_headers(self, cols):
        # Trim whitespace only; keep case for readability
        return [c.strip() if isinstance(c, str) else c for c in cols]

    def _detect_rel_type_col(self, file_path, requested_name: str | None):
        """
        Return the actual column name (case-insensitive) if present; else None.
        """
        try:
            hdr = pd.read_csv(file_path, nrows=0, dtype=str, keep_default_na=True)
        except Exception as e:
            self.logger.error("Failed reading headers for %s: %s", file_path, e)
            return None

        cols = self._normalize_headers(list(hdr.columns))
        # If the caller asked for a specific name, try that first (case-insensitive)
        candidates = []
        if requested_name:
            candidates.append(requested_name)
        # Also try some common aliases automatically
        candidates += ["activity", "relationship_type", "relation", "rel_type"]

        lower_map = {c.lower(): c for c in cols}
        for cand in candidates:
            lc = cand.strip().lower()
            if lc in lower_map:
                return lower_map[lc]
        return None


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
            escaped_value = (
                str(value)
                .replace("\\", "\\\\")
                .replace("'", "\\'")
                .replace('"', '\\"')
                .replace("\n", "\\\n")
            )
            return f"'{escaped_value}'"

    def _process_properties(
        self, row, source_column, destination_column, rel_type_column, standard_id
    ):
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
        properties = {f"`{k}`": v for k, v in properties.items()}  ##last_edit
        return properties

    def _generate_query(
        self,
        source,
        target,
        relationship_type,
        properties,
        source_label,
        destination_label,
        standard_id,
        source_column,
        destination_column,
    ):
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
        source_key = (
            standard_id[source_column]
            if source_column in standard_id
            else source_column
        )
        destination_key = (
            standard_id[destination_column]
            if destination_column in standard_id
            else destination_column
        )

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

    def generate_cypher_queries_from_file(
        self, file_path, rel_type, source_label, destination_label, rel_type_column=None
    ):
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
            "GeneID": "GeneID", "Gene ID": "GeneID", "Target GeneID": "GeneID",
            "geneids": "GeneID", "AssayID": "AssayID", "Assay ID": "AssayID",
            "AID": "AssayID", "aid": "AssayID", "ID_1": "CompoundID", "ID_2": "CompoundID",
            "substratecid": "CompoundID", "metabolitecid": "CompoundID",
            "Compound ID": "CompoundID", "CompoundID": "CompoundID", "CID": "CompoundID",
            "Similar CIDs": "CompoundID", "Target Accession": "ProteinRefSeqAccession",
            "geneid": "GeneID", "target_geneid": "GeneID", "cid": "CompoundID",
        }
        self.logger.info(f"Reading data from CSV file: {file_path}")

        # Read as strings to avoid dtype surprises; treat "__nan__" as NA
        df = pd.read_csv(
            file_path, dtype=str, keep_default_na=True, na_values=["__nan__"], low_memory=False
        )
        df.columns = self._normalize_headers(df.columns.tolist())
        df = df.dropna(axis=1, how="all")

        if df.shape[1] < 2:
            self.logger.error("CSV %s has <2 usable columns after cleaning.", file_path)
            return

        if rel_type == "CO_OCCURS_IN_LITERATURE":
            df.rename(
                columns={df.columns[0]: list(ast.literal_eval(df.iloc[0, 0]).keys())[0]},
                inplace=True,
            )

        source_column, destination_column = df.columns[:2]
        df = df.dropna(subset=[source_column, destination_column], how="any")
        if df.empty:
            self.logger.error("CSV %s contains no valid rows after filtering.", file_path)
            return

        # Only drop NA on rel_type_column if present in this file
        if rel_type_column and rel_type_column in df.columns:
            before = len(df)
            df = df.dropna(subset=[rel_type_column])
            if df.empty:
                self.logger.warning(
                    "No rows remain in %s after dropping NA on '%s' (started %d).",
                    file_path, rel_type_column, before
                )
                return
        else:
            # If column was requested but not found, fall back to constant
            if rel_type_column:
                self.logger.warning(
                    "Requested rel_type_column '%s' not found in %s. "
                    "Falling back to constant rel_type '%s'.",
                    rel_type_column, file_path, rel_type
                )
            rel_type_column = None

        # (Optional) write a small preview for debugging
        try:
            df.head(50).to_csv("Data/test.csv", index=False)
            self.logger.info("Wrote preview to Data/test.csv (first 50 rows).")
        except Exception as e:
            self.logger.debug("Could not write Data/test.csv: %s", e)

        for _, row in df.iterrows():
            src_raw = row[source_column]
            dst_raw = row[destination_column]
            relationship_type = (row[rel_type_column] if rel_type_column else rel_type) or rel_type
            relationship_type = str(relationship_type).replace("/", "OR")

            properties = self._process_properties(
                row, source_column, destination_column, rel_type_column, standard_id
            )

            if rel_type == "IS_SIMILAR_TO":
                targets = ast.literal_eval(dst_raw)
                for target in targets:
                    yield self._generate_query(
                        src_raw, target, relationship_type, properties,
                        source_label, destination_label, standard_id,
                        source_column, destination_column
                    )

            elif rel_type == "CO_OCCURS_IN_LITERATURE":
                src = ast.literal_eval(src_raw)
                if isinstance(src, dict):
                    src = list(src.values())[0]
                targets = ast.literal_eval(dst_raw)
                if isinstance(targets, dict):
                    for target in targets.values():
                        yield self._generate_query(
                            src, target, relationship_type, properties,
                            source_label, destination_label, standard_id,
                            source_column, destination_column
                        )

            elif rel_type == "ENCODES":
                yield self._generate_query(
                    str(src_raw), str(dst_raw), relationship_type, properties,
                    source_label, destination_label, standard_id,
                    source_column, destination_column
                )
            else:
                # Keep as strings (safer key matching in Neo4j)
                yield self._generate_query(
                    str(src_raw), str(dst_raw), relationship_type, properties,
                    source_label, destination_label, standard_id,
                    source_column, destination_column
                )

        self.logger.info("Cypher queries generated successfully.")

    def generate_cypher_queries_from_directories(
        self, directory, rel_type, source_label, destination_label, rel_type_column=None
    ):
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

            # Detect per-file relationship-type column (case-insensitive)
            per_file_rel_col = self._detect_rel_type_col(csv_file, rel_type_column)

            for query in self.generate_cypher_queries_from_file(
                csv_file, rel_type, source_label, destination_label, per_file_rel_col
            ):
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
                batch = queries[i : i + batch_size]
                try:
                    for query in batch:
                        session.run(query)
                        self.logger.debug(f"Executed query: {query}")
                    self.logger.info(
                        f"Executed batch {i // batch_size + 1}"
                        f"of {len(queries) // batch_size + 1}"
                    )
                except Exception as e:
                    self.logger.error(
                        "Failed to execute batch starting at query %s: %s", i, str(e)
                    )

        self.logger.info("All queries executed.")

    def process_and_add_relationships(
        self, file_path, rel_type, source_label, destination_label, rel_type_column=None
    ):
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
        self.logger.info("Processing and adding relationships from file: %s", file_path)
        queries = list(
            self.generate_cypher_queries_from_file(
                file_path, rel_type, source_label, destination_label, rel_type_column
            )
        )
        self.execute_queries(queries)
        self.logger.info(
            "Successfully processed and added relationships from file: %s", file_path
        )

    def process_and_add_relationships_from_directory(
        self,
        directory_path,
        rel_type,
        source_label,
        destination_label,
        rel_type_column=None,
    ):
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
        self.logger.info(
            "Processing and adding relationships from directory: %s", directory_path
        )
        queries = list(
            self.generate_cypher_queries_from_directories(
                directory_path,
                rel_type,
                source_label,
                destination_label,
                rel_type_column,
            )
        )
        if not queries:
            self.logger.error(
                "No valid relationships found in the directory %s.", directory_path
            )
            return

        self.execute_queries(queries)
        self.logger.info(
            "Successfully processed and added relationships from directory: %s",
            directory_path,
        )

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
                self.logger.info(
                    f"Bidirectional relationships for {relationship_type} created successfully."
                )
            except Exception as e:
                self.logger.error(
                    f"Failed to create bidirectional relationships for {relationship_type}: {str(e)}"
                )
