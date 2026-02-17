"""
Module for adding node data from CSV files to a Neo4j database.

Upgrades vs previous version:
- Uses parameterized Cypher + UNWIND batching (much faster than per-row string queries).
- Converts pandas/NumPy NaN to None (Neo4j cannot store NaN).
- Optionally normalizes multi-valued fields (e.g., *_IDs, Patent_IDs, DBURLs) into list properties.
- Avoids overwriting existing properties with nulls by default (skips None values).
- Keeps backward-compatible public API methods where possible.
"""

from __future__ import annotations

import logging
import math
import os
import re
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd

from chemgraphbuilder.neo4jdriver import Neo4jBase

logger = logging.getLogger(__name__)


# Heuristics for "multi-valued" columns (pipe-separated) that we should store as list properties in Neo4j
_PIPE_LIST_COL_RE = re.compile(r".*_IDs$|.*_ID(s)?$")  # conservative; we additionally gate by value containing '|'
_PIPE_LIST_COL_EXACT = {
    "DBURLs",
    "Patent_IDs",
    "MeSH_Entry_Terms",
    "MeSH_Pharmacological_Classification",
    "MeSH_Embedded_IDs",
    "MeSH_All_IDs",
    "EC_Numbers",
}
# Some fields are semicolon-separated lists
_SEMI_LIST_COL_EXACT = {"OntologySources"}


def _is_nan(x: Any) -> bool:
    try:
        return x is None or (isinstance(x, float) and math.isnan(x))
    except Exception:
        return False


def _to_native(x: Any) -> Any:
    """Convert pandas/NumPy scalars to Python scalars."""
    # Pandas uses numpy dtypes; converting via item() when available
    try:
        if hasattr(x, "item"):
            return x.item()
    except Exception:
        pass
    return x


class AddGraphNodes(Neo4jBase):
    """
    Add node data from CSV files to Neo4j with batching and type coercion.
    """

    def __init__(self, driver):
        super().__init__()
        self.driver = driver
        self.logger.info("AddGraphNodes initialized.")

    @staticmethod
    def create_uniqueness_constraint(driver, label: str, unique_property: str) -> None:
        query = (
            f"CREATE CONSTRAINT IF NOT EXISTS FOR (n:{label}) "
            f"REQUIRE n.{unique_property} IS UNIQUE"
        )
        with driver.session() as session:
            session.run(query)
        logging.info("Uniqueness constraint ensured: (%s.%s) UNIQUE", label, unique_property)

    @staticmethod
    def create_index(driver, label: str, prop: str) -> None:
        query = f"CREATE INDEX IF NOT EXISTS FOR (n:{label}) ON (n.{prop})"
        with driver.session() as session:
            session.run(query)
        logging.info("Index ensured: (%s.%s)", label, prop)

    def _coerce_value(
        self,
        value: Any,
        col_name: str,
        *,
        normalize_lists: bool = True,
    ) -> Any:
        """
        Convert values to Neo4j-compatible property types:
        - None instead of NaN
        - Python scalars (int/float/bool/str)
        - list[str] for multi-valued columns (pipe/semicolon separated), if enabled
        """
        value = _to_native(value)

        if _is_nan(value):
            return None

        # Keep bool as bool
        if isinstance(value, bool):
            return value

        # Convert ints/floats
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            # Neo4j can't store NaN; already handled above
            if float(value).is_integer():
                return int(value)
            return float(value)

        # Normalize strings
        if isinstance(value, str):
            s = value.strip()
            if s == "":
                return ""  # keep empty string as-is

            if normalize_lists:
                # Semicolon list (OntologySources)
                if col_name in _SEMI_LIST_COL_EXACT and ";" in s:
                    parts = [p.strip() for p in s.split(";") if p.strip()]
                    return parts if len(parts) > 1 else parts[0]

                # Pipe-separated lists for known/likely multi-valued columns
                if ("|" in s) and (_PIPE_LIST_COL_RE.match(col_name) or col_name in _PIPE_LIST_COL_EXACT):
                    parts = [p.strip() for p in s.split("|") if p.strip()]
                    return parts if len(parts) > 1 else parts[0]

            return s

        # Fallback: stringify unknown types
        try:
            return str(value)
        except Exception:
            return None

    def _coerce_row(
        self,
        row: Dict[str, Any],
        *,
        skip_nulls: bool = True,
        normalize_lists: bool = True,
    ) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        for k, v in row.items():
            vv = self._coerce_value(v, k, normalize_lists=normalize_lists)
            if skip_nulls and vv is None:
                continue
            out[k] = vv
        return out

    def read_csv_file(self, file_path: str) -> pd.DataFrame:
        self.logger.info("Reading CSV: %s", file_path)
        return pd.read_csv(file_path)

    def upsert_nodes_from_dataframe(
        self,
        df: pd.DataFrame,
        *,
        label: str,
        unique_property: str,
        batch_size: int = 50,
        skip_nulls: bool = True,
        normalize_lists: bool = True,
        create_unique_constraint: bool = True,
        create_unique_index: bool = True,
    ) -> int:
        """
        Upsert nodes into Neo4j from a dataframe using UNWIND batches.

        Returns number of rows attempted (after dropping missing unique_property).
        """
        if unique_property not in df.columns:
            raise ValueError(f"Unique property '{unique_property}' not in dataframe columns for label '{label}'")

        total = len(df)
        df = df.dropna(subset=[unique_property], how="any")
        dropped = total - len(df)
        if dropped:
            self.logger.warning("Dropped %d rows with missing %s for label %s", dropped, unique_property, label)

        if create_unique_constraint:
            self.create_uniqueness_constraint(self.driver, label, unique_property)
        if create_unique_index:
            # A separate index is usually redundant when a uniqueness constraint exists,
            # but we keep it here for compatibility with older Neo4j setups.
            try:
                self.create_index(self.driver, label, unique_property)
            except Exception:
                pass

        # Build batches
        rows_attempted = len(df)
        self.logger.info("Upserting %d %s nodes (batch_size=%d)...", rows_attempted, label, batch_size)

        query = f"""
        UNWIND $rows AS row
        MERGE (n:{label} {{ {unique_property}: row.{unique_property} }})
        SET n += row
        """

        # Neo4j transactions should not be too big; use batches
        with self.driver.session() as session:
            for start in range(0, rows_attempted, batch_size):
                chunk = df.iloc[start : start + batch_size]
                rows: List[Dict[str, Any]] = []
                for _, r in chunk.iterrows():
                    row_dict = self._coerce_row(
                        r.to_dict(),
                        skip_nulls=skip_nulls,
                        normalize_lists=normalize_lists,
                    )
                    # Ensure unique_property exists (and is not None due to coercion)
                    if unique_property in row_dict and row_dict[unique_property] is not None:
                        rows.append(row_dict)

                if not rows:
                    continue

                session.execute_write(lambda tx: tx.run(query, rows=rows))

        self.logger.info("Upsert complete for label %s (%d rows).", label, rows_attempted)
        return rows_attempted

    def upsert_nodes_from_csv(
        self,
        file_path: str,
        *,
        label: str,
        unique_property: str,
        batch_size: int = 50,
        skip_nulls: bool = True,
        normalize_lists: bool = True,
    ) -> int:
        df = self.read_csv_file(file_path)
        return self.upsert_nodes_from_dataframe(
            df,
            label=label,
            unique_property=unique_property,
            batch_size=batch_size,
            skip_nulls=skip_nulls,
            normalize_lists=normalize_lists,
        )

    # -------------------------------------------------------------------------
    # Backward-compatible wrapper methods (kept to minimize breakage)
    # -------------------------------------------------------------------------
    def process_and_add_nodes(self, file_path: str, label: str, unique_property: str) -> None:
        """Legacy API: now uses fast upsert with defaults."""
        self.upsert_nodes_from_csv(
            file_path,
            label=label,
            unique_property=unique_property,
            batch_size=50,
            skip_nulls=True,
            normalize_lists=True,
        )

    def process_and_add_nodes_from_directory(self, directory_path: str, label: str, unique_property: str) -> None:
        """Legacy API: concatenate CSVs in directory then upsert."""
        dfs = [
            pd.read_csv(os.path.join(directory_path, f))
            for f in os.listdir(directory_path)
            if f.endswith(".csv")
        ]
        if not dfs:
            self.logger.warning("No CSVs found in directory: %s", directory_path)
            return
        combined = pd.concat(dfs, ignore_index=True)
        self.upsert_nodes_from_dataframe(
            combined,
            label=label,
            unique_property=unique_property,
            batch_size=50,
            skip_nulls=True,
            normalize_lists=True,
        )
