"""
Module to load node data into a Neo4j graph database for different node types.

Upgrades vs previous version:
- Supports all node types needed by the thesis schema:
  Compound, BioAssay, Protein, Gene, ExperimentalContext, AssayEndpoint.
- Prefers *_WithOntologies.csv when available (keeps processed as fallback).
- Uses fast UNWIND batching via AddGraphNodes.upsert_nodes_from_csv().
- Cleans NaN/None correctly and can normalize list-like columns into Neo4j list properties.
- Can load multiple labels in one run.
- Works even when config/node_schema.yml is missing (uses safe defaults).

This script loads NODES only. Relationships can be loaded separately from AllDataConnected.csv
(or you can extend this script with a relationship loader).
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Dict, List, Optional

import pandas as pd
import yaml
from neo4j import GraphDatabase

from chemgraphbuilder.add_graph_nodes import AddGraphNodes

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Fallback mapping if YAML schema is missing or incomplete.
# "label" is the Neo4j label; "file" is the default CSV filename (prefer WithOntologies).
DEFAULT_LABEL_MAPPING: Dict[str, Dict[str, str]] = {
    "Compound": {
        "unique_property": "CompoundID",
        "file": "Compound_Properties_WithOntologies.csv",
    },
    "BioAssay": {
        "unique_property": "AssayID",
        "file": "Assay_Properties_WithOntologies.csv",
    },
    # allow a common alias
    "Assay": {
        "unique_property": "AssayID",
        "file": "Assay_Properties_WithOntologies.csv",
    },
    "Gene": {
        "unique_property": "GeneID",
        "file": "Gene_Properties_WithOntologies.csv",
    },
    "Protein": {
        "unique_property": "ProteinRefSeqAccession",
        "file": "Protein_Properties_WithOntologies.csv",
    },
    "ExperimentalContext": {
        "unique_property": "ExperimentalContextID",
        "file": "ExperimentalContext_Properties_WithOntologies.csv",
    },
    "AssayEndpoint": {
        "unique_property": "AssayEndpointID",
        "file": "AssayEndpoint_Properties_WithOntologies.csv",
    },
}


class GraphNodesLoader:
    """
    Loads CSV node tables into Neo4j according to a schema mapping (YAML) or fallback defaults.
    """

    def __init__(
        self,
        uri: str,
        username: str,
        password: str,
        *,
        schema_path: str = "config/node_schema.yml",
        data_dir: str = "Data/Nodes",
    ):
        self.driver = GraphDatabase.driver(uri, auth=(username, password))
        self.node_data_adder = AddGraphNodes(self.driver)
        self.data_dir = data_dir

        self.label_mapping = self._load_schema_mapping(schema_path)
        logger.info("GraphNodesLoader initialized (data_dir=%s).", self.data_dir)

    def close(self) -> None:
        try:
            self.driver.close()
        except Exception:
            pass

    def _resolve_path(self, p: str) -> str:
        if os.path.isabs(p):
            return p
        # If schema provides a relative path, resolve relative to CWD.
        # If the file is not found, try resolving under data_dir.
        if os.path.exists(p):
            return p
        candidate = os.path.join(self.data_dir, os.path.basename(p))
        return candidate

    def _load_schema_mapping(self, schema_path: str) -> Dict[str, Dict[str, str]]:
        mapping: Dict[str, Dict[str, str]] = {}

        if schema_path and os.path.exists(schema_path):
            try:
                with open(schema_path, "r", encoding="utf-8") as f:
                    schema = yaml.safe_load(f) or {}
                nodes = (schema.get("nodes") or {})
                for _, node_cfg in nodes.items():
                    label = node_cfg.get("label")
                    unique_property = node_cfg.get("unique_property")
                    # prefer ontologies file if present in schema
                    fp = node_cfg.get("csv_ontologies") or node_cfg.get("csv_with_ontologies") or node_cfg.get("csv_processed")
                    if label and unique_property and fp:
                        mapping[label] = {
                            "unique_property": unique_property,
                            "file_path": self._resolve_path(fp),
                        }
            except Exception as exc:
                logger.warning("Failed to load schema yaml (%s): %s. Falling back to defaults.", schema_path, exc)

        # Merge fallback defaults (do not override schema-provided entries)
        for label, cfg in DEFAULT_LABEL_MAPPING.items():
            if label not in mapping:
                fp = os.path.join(self.data_dir, cfg["file"])
                mapping[label] = {"unique_property": cfg["unique_property"], "file_path": fp}

        return mapping

    def validate_node_type(self, label: str) -> bool:
        """
        Validate CSV for a given label:
        - exists
        - contains unique_property
        - report duplicates
        """
        if label not in self.label_mapping:
            logger.error("No mapping found for label: %s", label)
            return False

        unique_property = self.label_mapping[label]["unique_property"]
        file_path = self.label_mapping[label]["file_path"]

        file_path = self._resolve_path(file_path)
        if not os.path.exists(file_path):
            logger.error("CSV file for label %s not found: %s", label, file_path)
            return False

        try:
            df = pd.read_csv(file_path)
        except Exception as exc:
            logger.error("Failed to read CSV for label %s (%s): %s", label, file_path, exc)
            return False

        if unique_property not in df.columns:
            logger.error("Unique property '%s' not in %s for label %s", unique_property, file_path, label)
            return False

        dup_count = int(df[unique_property].duplicated().sum())
        if dup_count > 0:
            logger.warning("Found %d duplicate %s values in %s", dup_count, unique_property, file_path)

        logger.info("Validation OK for label %s (file: %s)", label, file_path)
        return True

    def load_data_for_node_type(
        self,
        label: str,
        *,
        batch_size: int = 50,
        normalize_lists: bool = True,
        skip_nulls: bool = True,
    ) -> None:
        if label not in self.label_mapping:
            logger.error("No mapping found for label: %s", label)
            return

        unique_property = self.label_mapping[label]["unique_property"]
        file_path = self._resolve_path(self.label_mapping[label]["file_path"])

        if not os.path.exists(file_path):
            logger.error("File not found for label %s: %s", label, file_path)
            return

        # Ensure uniqueness constraint
        self.node_data_adder.create_uniqueness_constraint(self.driver, label=label, unique_property=unique_property)

        # Upsert nodes
        self.node_data_adder.upsert_nodes_from_csv(
            file_path,
            label=label,
            unique_property=unique_property,
            batch_size=batch_size,
            normalize_lists=normalize_lists,
            skip_nulls=skip_nulls,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Load node tables into Neo4j (WithOntologies preferred).")
    parser.add_argument("--uri", required=True, help="URI for the Neo4j database (e.g. bolt://localhost:7687)")
    parser.add_argument("--username", required=True, help="Neo4j username")
    parser.add_argument("--password", required=True, help="Neo4j password")

    # Backwards-compatible: --label (single) or --labels (multiple)
    parser.add_argument("--label", default=None, help="Single label to load (legacy).")
    parser.add_argument(
        "--labels",
        default=None,
        help="Comma-separated labels to load (recommended). If omitted, loads all known labels.",
    )

    parser.add_argument(
        "--schema_path",
        default="config/node_schema.yml",
        help="Path to node schema YAML (optional). Default: config/node_schema.yml",
    )
    parser.add_argument(
        "--data_dir",
        default=os.path.join("Data", "Nodes"),
        help="Base directory for node CSVs (used for defaults + resolving missing files). Default: Data/Nodes",
    )
    parser.add_argument("--batch_size", type=int, default=50, help="Batch size for UNWIND upserts. Default: 500")
    parser.add_argument(
        "--no_normalize_lists",
        action="store_true",
        help="Disable splitting pipe/semicolon-separated columns into Neo4j list properties.",
    )
    parser.add_argument(
        "--overwrite_with_nulls",
        action="store_true",
        help="If set, include null values in SET n += row (may remove existing properties). Default: skip nulls.",
    )

    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Validate CSV(s) only; do not write to Neo4j.",
    )

    args = parser.parse_args()

    loader = GraphNodesLoader(
        args.uri,
        args.username,
        args.password,
        schema_path=args.schema_path,
        data_dir=args.data_dir,
    )

    try:
        if args.label and args.labels:
            logger.warning("Both --label and --labels provided; using --labels.")

        if args.labels:
            labels = [x.strip() for x in args.labels.split(",") if x.strip()]
        elif args.label:
            labels = [args.label.strip()]
        else:
            # load all labels defined in mapping (defaults + schema)
            labels = list(loader.label_mapping.keys())

        normalize_lists = not args.no_normalize_lists
        skip_nulls = not args.overwrite_with_nulls

        if args.dry_run:
            ok = True
            for lb in labels:
                ok = loader.validate_node_type(lb) and ok
            if ok:
                logger.info("✅ Dry-run validation passed for labels: %s", labels)
            else:
                logger.error("❌ Dry-run validation failed for at least one label.")
            return

        for lb in labels:
            if not loader.validate_node_type(lb):
                logger.error("Skipping label %s due to validation failure.", lb)
                continue
            loader.load_data_for_node_type(
                lb,
                batch_size=args.batch_size,
                normalize_lists=normalize_lists,
                skip_nulls=skip_nulls,
            )

        logger.info("✅ Node loading finished for labels: %s", labels)

    finally:
        loader.close()


if __name__ == "__main__":
    main()
