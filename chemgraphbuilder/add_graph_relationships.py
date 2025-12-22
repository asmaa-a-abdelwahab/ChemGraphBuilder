"""
Module for adding relationship data from CSV files to a Neo4j database.

This module provides a class and methods to read relationship data from CSV files
and add them to a Neo4j database, including generating Cypher queries.
"""

import glob
import ast
import re
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
    # @staticmethod
    # def _sanitize_rel_type(name: str) -> str:
    #     """
    #     Sanitizes a relationship type name by converting it to uppercase, replacing
    #     non-alphanumeric characters with underscores, and stripping trailing underscores.
    #     If the resulting string does not start with a letter or underscore, it
    #     is prefixed with an underscore.
    #     """
        
    #     s = str(name).upper()
    #     s = re.sub(r"[^A-Z0-9_]", "_", s)
    #     s = re.sub(r"_+", "_", s).strip("_")
    #     return s

    @staticmethod
    def _pick_id_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
        cols_lower = {c.lower(): c for c in df.columns}
        for cand in candidates:
            if cand.lower() in cols_lower:
                return cols_lower[cand.lower()]
        return None

    @staticmethod
    def _id_candidates_for_label(label: str) -> list[str]:
        label = (label or "").lower()
        if label in ("bioassay", "assay"):
            return ["AssayID", "AID", "aid"]
        if label in ("compound","Compound"):
            return ["CompoundID", "CID", "cid"]
        if label in ("gene","Gene"):
            return ["GeneID", "geneid", "Target GeneID", "target_geneid", "GeneSymbol", "genesymbol"]
        if label in ("protein",):
            return ["ProteinRefSeqAccession", "Target Accession", "target_accession"]
        return ["AssayID", "AID", "CID", "GeneID", "geneid"]


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

        # ✅ Use two MATCH clauses → no Cartesian product warning
        query = (
            f"MATCH (a:{source_label} {{{source_key}: {source_value}}})\n"
            f"MATCH (b:{destination_label} {{{destination_key}: {target_value}}})\n"
            f"MERGE (a)-[r:{relationship_type}]->(b)"
        )

        if properties:
            set_clauses = [
                f"r.{prop.replace(' ', '').replace('[', '_').replace(']', '_')}="
                f"{self._generate_property_string(value)}"
                for prop, value in properties.items()
            ]
            if set_clauses:
                query += " SET " + ", ".join(set_clauses)

        self.logger.debug("Generated query: %s", query)
        return query

    def _id_key_for_label(self, label: str) -> str:
        """
        Map a node label to the key we expect inside the ID_1/ID_2 dicts.
        Only keys used by your co-occurrence files are handled here.
        """
        l = (label or "").strip().lower()
        if l in ("compound",):
            return "CID"
        if l in ("gene",):
            return "GeneSymbol"
        if l in ("bioassay", "assay"):
            return "AID"
        if l in ("protein",):
            return "Target Accession"
        # default fallback; unlikely used for co-occurrence
        return "CID"


    def _normalize_cooccurrence(self, df: pd.DataFrame, source_label: str, destination_label: str) -> tuple[pd.DataFrame, str, str]:
        """
        Normalize cpd–cpd and cpd–gene co-occurrence rows into a common, 2-column form:
        <src_col>, <dst_col>, plus optional 'Evidence' if present.
        - Assumes the raw file has columns ID_1, ID_2 (each a dict-like string).
        - Picks values by the expected keys inferred from node labels.
        Returns: (normalized_df, src_col_name, dst_col_name)
        """
        if "ID_1" not in df.columns or "ID_2" not in df.columns:
            src_key = self._id_key_for_label(source_label)
            dst_key = self._id_key_for_label(destination_label)
            return pd.DataFrame(columns=[src_key + "_src", dst_key + "_dst", "Evidence"]), src_key + "_src", dst_key + "_dst"

        def to_dict(x):
            try:
                return ast.literal_eval(x) if isinstance(x, str) else {}
            except Exception:
                return {}

        src_key = self._id_key_for_label(source_label)      # "CID" or "GeneSymbol"
        dst_key = self._id_key_for_label(destination_label) # "CID" or "GeneSymbol"

        d1 = df["ID_1"].apply(to_dict)
        d2 = df["ID_2"].apply(to_dict)

        src_vals = d1.apply(lambda d: d.get(src_key)).fillna(d2.apply(lambda d: d.get(src_key)))
        dst_vals = d2.apply(lambda d: d.get(dst_key)).fillna(d1.apply(lambda d: d.get(dst_key)))

        src_col = f"{src_key}_src"
        dst_col = f"{dst_key}_dst"

        out = pd.DataFrame({
            src_col: src_vals.astype(str),
            dst_col: dst_vals.astype(str),
        })

        if "Evidence" in df.columns:
            out["Evidence"] = df["Evidence"]

        out = out.dropna(subset=[src_col, dst_col], how="any")

        return out, src_col, dst_col


    def _parse_listish(self, value):
        """
        Return a list of string targets from a variety of formats:
        - "[123, 456]" or "['123','456']"  -> ["123","456"]
        - "123,456" or "123 ; 456"         -> ["123","456"]
        - "12345" (scalar)                 -> ["12345"]
        - NaN / empty                      -> []
        """
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return []

        # Work with strings thereafter
        s = str(value).strip()
        if not s or s.lower() in {"nan", "none"}:
            return []

        # Looks like a Python/JSON list?
        if (s.startswith("[") and s.endswith("]")) or (s.startswith("(") and s.endswith(")")):
            try:
                v = ast.literal_eval(s)
                if isinstance(v, (list, tuple, set)):
                    return [str(x).strip() for x in v if str(x).strip()]
                # If it’s a scalar inside the brackets for some reason
                return [str(v)]
            except Exception:
                pass  # fall through to split heuristics

        # Comma/semicolon separated?
        if ("," in s) or (";" in s) or ("|" in s):
            # Normalize multiple separators
            for sep in [";", "|"]:
                s = s.replace(sep, ",")
            return [t.strip() for t in s.split(",") if t.strip()]

        # Otherwise treat as single value
        return [s]


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
            "geneid": "GeneID", "target_geneid": "GeneID", "cid": "CompoundID",     "GeneSymbol": "GeneSymbol",
            "genesymbol": "GeneSymbol", "CID": "CompoundID", "cid": "CompoundID",
            "LinkID": "CompoundID", "linkid": "CompoundID", "ID": "CompoundID",
            "CID_src": "CompoundID", "CID_dst": "CompoundID", "GeneSymbol_src": "GeneSymbol",
            "GeneSymbol_dst": "GeneSymbol", "AID_src": "AssayID", "AID_dst": "AssayID",
        }
        # self.logger.info(f"Reading data from CSV file: {file_path}")

        # Read as strings to avoid dtype surprises; treat "__nan__" as NA
        df = pd.read_csv(
            file_path,
            dtype=str,
            keep_default_na=True,
            na_values=["__nan__"],
            engine="python",           # tolerant to weird quoting/newlines
            on_bad_lines="skip",       # skip broken Evidence rows
            quotechar='"',
            escapechar='\\',
            index_col=False
        )

        df.columns = self._normalize_headers(df.columns.tolist())
        df = df.dropna(axis=1, how="all")

        if df.shape[1] < 2:
            self.logger.error("CSV %s has <2 usable columns after cleaning.", file_path)
            return

        # --- 1) Co-occurrence (Cpd–Cpd, Cpd–Gene) ---
        if rel_type in {
            "CO_OCCURS_IN_LITERATURE",
            "Compound_Gene_CoOccurrence",
            "Compound_Compound_CoOccurrence",
        }:
            # Normalize (works for both CID↔CID and GeneSymbol↔CID)
            df, source_column, destination_column = self._normalize_cooccurrence(
                df, source_label, destination_label
            )
            if df.empty:
                self.logger.error(
                    "CSV %s contains no valid co-occurrence rows after normalization.",
                    file_path,
                )
                return

        # --- 2) Compound transformation (substrate → metabolite) ---
        elif (
            # Normal case: explicit relationship type
            rel_type in {"Compound_Transformation", "METABOLIZED_TO", "TRANSFORMS_TO"}
            # Fallback: detect by columns, in case rel_type was sanitized
            or (
                source_label.lower() == "compound"
                and destination_label.lower() == "compound"
                and {"substratecid", "metabolitecid"}.issubset(df.columns)
            )
        ):
            if {"substratecid", "metabolitecid"}.issubset(df.columns):
                source_column = "substratecid"
                destination_column = "metabolitecid"
            else:
                self.logger.error(
                    "Compound_Transformation loader: expected 'substratecid' "
                    "and 'metabolitecid' in %s",
                    file_path,
                )
                return

        # --- 3) Compound similarity (CID → [similar CIDs]) ---
        elif rel_type in {"IS_SIMILAR_TO", "Compound_Similarity", "COMPOUND_SIMILARITY"}:
            list_col_candidates = [
                "Similar CIDs",
                "Similar_CIDs",
                "similar_cids",
                "similarCIDs",
                "similar_ids",
            ]
            src_candidates = ["CID", "cid", "CompoundID", "Compound ID", "LinkID", "ID"]
            source_column = self._pick_id_column(df, src_candidates)
            destination_column = next(
                (c for c in list_col_candidates if c in df.columns), None
            )

            if not source_column or not destination_column:
                self.logger.error(
                    "Similarity loader: could not find required columns in %s "
                    "(source one of %s, destination one of %s).",
                    file_path,
                    src_candidates,
                    list_col_candidates,
                )
                return

            df = df.dropna(subset=[source_column], how="any")
            if df.empty:
                self.logger.error(
                    "CSV %s contains no valid rows after filtering (source).",
                    file_path,
                )
                return

        # --- 4) Generic relationships (everything else) ---
        else:
            src_candidates = self._id_candidates_for_label(source_label)
            dst_candidates = self._id_candidates_for_label(destination_label)

            source_column = self._pick_id_column(df, src_candidates)
            destination_column = self._pick_id_column(df, dst_candidates)

            if not source_column or not destination_column:
                self.logger.error(
                    "Could not detect ID columns for %s → %s in %s",
                    source_label,
                    destination_label,
                    file_path,
                )
                return

            df = df.dropna(subset=[source_column, destination_column], how="any")
            if df.empty:
                self.logger.error(
                    "CSV %s contains no valid rows after filtering.", file_path
                )
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

        # # (Optional) write a small preview for debugging
        # try:
        #     df.head(50).to_csv("Data/test.csv", index=False)
        #     self.logger.info("Wrote preview to Data/test.csv (first 50 rows).")
        # except Exception as e:
        #     self.logger.debug("Could not write Data/test.csv: %s", e)

        for _, row in df.iterrows():
            src_raw = row[source_column]
            dst_raw = row[destination_column]
            relationship_type = (row[rel_type_column] if rel_type_column else rel_type) or rel_type
            relationship_type = str(relationship_type).replace("/", "OR")


            properties = self._process_properties(
                row, source_column, destination_column, rel_type_column, standard_id
            )

            if rel_type == "IS_SIMILAR_TO":
                targets = self._parse_listish(row.get(destination_column))
                if not targets:
                    continue
                for target in targets:
                    # skip degenerate self-links early (we also have a cleanup method)
                    if str(target) == str(src_raw):
                        continue
                    yield self._generate_query(
                        str(src_raw), str(target), relationship_type, properties,
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

            # Only auto-detect when rel_type is None (dynamic relationship types)
            if rel_type is None:
                per_file_rel_col = self._detect_rel_type_col(csv_file, rel_type_column)
            else:
                per_file_rel_col = rel_type_column   # keep whatever was passed (None)

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
