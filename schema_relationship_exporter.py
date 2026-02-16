from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


# -----------------------------
# Helpers
# -----------------------------
def _ensure_dir(p: str | Path) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _safe_int(x: Any) -> int | None:
    if x is None:
        return None
    try:
        if pd.isna(x):
            return None
    except Exception:
        pass
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def _safe_str(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def _parse_py_dict(x: Any) -> dict:
    """
    Your LinkDB files store ID_1 / ID_2 / Evidence as python-dict-like strings,
    e.g. "{'CID': 1684}" or "{'ChemicalNeighbor': {...}}".
    """
    s = _safe_str(x)
    if not s:
        return {}
    try:
        d = ast.literal_eval(s)
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


def _split_pmids(x: Any) -> list[str]:
    s = _safe_str(x)
    if not s:
        return []
    # pmids may be "123|456" or "123;456"
    return sorted({m for m in re.findall(r"\d+", s)})


def _split_dois(x: Any) -> list[str]:
    s = _safe_str(x)
    if not s:
        return []
    # your Compound_Transformation.csv uses "|" delimiter
    parts = [p.strip() for p in s.split("|")]
    return sorted({p for p in parts if p})


def _join_list(values: Iterable[str], sep: str = ";", limit: int = 20000) -> str:
    out = sep.join(values)
    return out[:limit]


def _union_semicolon_lists(series: pd.Series) -> str:
    s = set()
    for v in series.dropna().astype(str):
        v = v.strip()
        if not v:
            continue
        for part in v.split(";"):
            part = part.strip()
            if part:
                s.add(part)
    return _join_list(sorted(s))


def _classify_inhibits(activity_name: Any) -> bool:
    """
    Based on your actual Compound_Gene_Relationship.csv values:
    IC50, Inhibition, INH => INHIBITS
    Km => INTERACTS_WITH (mode/kinetics)
    """
    s = _safe_str(activity_name).lower()
    return any(k in s for k in ["ic50", "inhib", "inh", "ki", "pki", "pic50"])


def _confidence(outcome: Any, value_um: Any) -> float | None:
    o = _safe_str(outcome).lower()
    has_val = _safe_str(value_um) != ""
    if o == "active" and has_val:
        return 0.9
    if o == "active":
        return 0.7
    if o == "inactive":
        return 0.3
    return None


def _parse_linkdb_block(evidence_str: Any, key: str) -> dict[str, Any]:
    """
    Evidence example:
      {'ChemicalNeighbor': {'ArticleCount':..., 'CooccurrenceScore':..., 'Article':[{'PMID':..., 'DOI':...}, ...]}}
    We keep: cooccurrence_score, article_count, pmids, dois, source
    """
    ev = _parse_py_dict(evidence_str)
    block = ev.get(key, {})
    if not isinstance(block, dict):
        return {}

    art_count = block.get("ArticleCount")
    score = block.get("CooccurrenceScore")
    articles = block.get("Article", [])
    pmids: set[str] = set()
    dois: set[str] = set()

    if isinstance(articles, list):
        for a in articles:
            if not isinstance(a, dict):
                continue
            pmid = _safe_str(a.get("PMID"))
            if pmid and pmid.lower() != "nan":
                # PMID sometimes numeric
                try:
                    pmids.add(str(int(float(pmid))))
                except Exception:
                    pmids.add(pmid)
            doi = _safe_str(a.get("DOI"))
            if doi and doi.lower() != "nan":
                dois.add(doi)

    return {
        "cooccurrence_score": _safe_int(score),
        "article_count": _safe_int(art_count),
        "pmids": _join_list(sorted(pmids)),
        "dois": _join_list(sorted(dois)),
        "source": "PubChemLinkDB",
    }


@dataclass
class SchemaExportConfig:
    # Inputs
    all_data_connected: str
    compound_gene_relationship: str
    linkdb_compound_compound_files: list[str]  # e.g., CID_*.csv
    linkdb_gene_symbol_files: list[str]        # e.g., Cpd_Gene_CoOccurrence_*.csv
    sdq_consolidatedcompoundtarget_files: list[str]  # Compound_Gene_Interaction_Outside_PubChem_*.csv
    compound_transformation: str               # Compound_Transformation.csv

    # Output
    out_dir: str = "Data/Relationships/SCHEMA"


class SchemaRelationshipExporter:
    def __init__(self, cfg: SchemaExportConfig):
        self.cfg = cfg
        self.out_dir = _ensure_dir(cfg.out_dir)

    # -----------------------------
    # CORE (from AllDataConnected.csv)
    # -----------------------------
    def export_core_from_all_data(self) -> None:
        df = pd.read_csv(self.cfg.all_data_connected)

        # 1) Compound -[:HAS_RESULT_IN]-> BioAssay
        cols = [
            "CID", "AID",
            "Activity Outcome", "Activity Name", "Activity Value [uM]",
            "Assay Name", "Assay Type", "PubMed ID", "RNAi",
            "SourceDatabase", "SourceEndpoint", "RetrievedAt",
        ]
        d = df[cols].copy()
        d["CID"] = d["CID"].apply(_safe_int)
        d["AID"] = d["AID"].apply(_safe_int)
        d = d.dropna(subset=["CID", "AID"]).drop_duplicates()

        d = d.rename(columns={
            "Activity Outcome": "outcome",
            "Activity Name": "effect_type",
            "Activity Value [uM]": "value_uM",
            "Assay Name": "assay_name",
            "Assay Type": "assay_type",
            "PubMed ID": "pubmed_id",
            "SourceDatabase": "source_db",
            "SourceEndpoint": "source_endpoint",
            "RetrievedAt": "retrieved_at",
        })
        d.to_csv(self.out_dir / "Compound_HAS_RESULT_IN_BioAssay.csv", index=False)

        # 2) Compound -[:PERTURBAGEN_OF]-> BioAssay (membership)
        d2 = df[["CID", "AID", "SourceDatabase", "SourceEndpoint", "RetrievedAt"]].copy()
        d2["CID"] = d2["CID"].apply(_safe_int)
        d2["AID"] = d2["AID"].apply(_safe_int)
        d2 = d2.dropna(subset=["CID", "AID"]).drop_duplicates()
        d2 = d2.rename(columns={
            "SourceDatabase": "source_db",
            "SourceEndpoint": "source_endpoint",
            "RetrievedAt": "retrieved_at",
        })
        d2.to_csv(self.out_dir / "Compound_PERTURBAGEN_OF_BioAssay.csv", index=False)

        # 3) BioAssay -[:HAS_TARGET_GENE]-> Gene
        d3 = df[["AID", "Target GeneID", "Activity Name", "SourceDatabase", "RetrievedAt"]].copy()
        d3["AID"] = d3["AID"].apply(_safe_int)
        d3["GeneID"] = d3["Target GeneID"].apply(_safe_int)
        d3 = d3.dropna(subset=["AID", "GeneID"]).drop_duplicates()
        d3 = d3.rename(columns={
            "Activity Name": "activity_name",
            "SourceDatabase": "source_db",
            "RetrievedAt": "retrieved_at",
        })[["AID", "GeneID", "activity_name", "source_db", "retrieved_at"]]
        d3.to_csv(self.out_dir / "BioAssay_HAS_TARGET_GENE_Gene.csv", index=False)

        # 4) BioAssay -[:HAS_TARGET_PROTEIN]-> Protein
        d4 = df[["AID", "Target Accession", "Activity Name", "SourceDatabase", "RetrievedAt"]].copy()
        d4["AID"] = d4["AID"].apply(_safe_int)
        d4["ProteinAccession"] = d4["Target Accession"].astype(str).str.strip()
        d4 = d4.dropna(subset=["AID", "ProteinAccession"]).drop_duplicates()
        d4 = d4.rename(columns={
            "Activity Name": "activity_name",
            "SourceDatabase": "source_db",
            "RetrievedAt": "retrieved_at",
        })[["AID", "ProteinAccession", "activity_name", "source_db", "retrieved_at"]]
        d4.to_csv(self.out_dir / "BioAssay_HAS_TARGET_PROTEIN_Protein.csv", index=False)

        # 5) Gene -[:ENCODES]-> Protein
        d5 = df[["Target GeneID", "Target Accession", "SourceDatabase", "RetrievedAt"]].copy()
        d5["GeneID"] = d5["Target GeneID"].apply(_safe_int)
        d5["ProteinAccession"] = d5["Target Accession"].astype(str).str.strip()
        d5 = d5.dropna(subset=["GeneID", "ProteinAccession"]).drop_duplicates()
        d5 = d5.rename(columns={
            "SourceDatabase": "source_db",
            "RetrievedAt": "retrieved_at",
        })[["GeneID", "ProteinAccession", "source_db", "retrieved_at"]]
        d5.to_csv(self.out_dir / "Gene_ENCODES_Protein.csv", index=False)

    # -----------------------------
    # Compound -> Protein (INHIBITS / INTERACTS_WITH)
    # -----------------------------
    def export_compound_protein_edges(self) -> None:
        df = pd.read_csv(self.cfg.compound_gene_relationship)

        df["CID"] = df["CID"].apply(_safe_int)
        df["ProteinAccession"] = df["Target Accession"].astype(str).str.strip()
        df["confidence"] = df.apply(
            lambda r: _confidence(r.get("Activity Outcome"), r.get("Activity Value [uM]")),
            axis=1,
        )
        df["is_inhibits"] = df["Activity Name"].apply(_classify_inhibits)

        base = df[[
            "CID", "ProteinAccession",
            "Activity Outcome", "Activity Name", "Activity Value [uM]",
            "confidence"
        ]].dropna(subset=["CID", "ProteinAccession"]).drop_duplicates()

        # INHIBITS
        inhib = base[base["Activity Name"].apply(_classify_inhibits)].copy()
        inhib = inhib.rename(columns={
            "Activity Outcome": "outcome",
            "Activity Name": "effect_type",
            "Activity Value [uM]": "value_uM",
        })
        inhib["source"] = "PubChem"
        inhib.to_csv(self.out_dir / "Compound_INHIBITS_Protein.csv", index=False)

        # INTERACTS_WITH
        inter = base[~base["Activity Name"].apply(_classify_inhibits)].copy()
        inter = inter.rename(columns={
            "Activity Outcome": "outcome",
            "Activity Name": "mode",
            "Activity Value [uM]": "value_uM",
        })
        inter["source"] = "PubChem"
        inter.to_csv(self.out_dir / "Compound_INTERACTS_WITH_Protein.csv", index=False)

    # -----------------------------
    # NEW: COOCCURS_IN_LITERATURE
    #   (:Compound)-[:COOCCURS_IN_LITERATURE]->(:Compound)
    #   (:GeneSymbol)-[:COOCCURS_IN_LITERATURE]->(:Compound)
    # -----------------------------
    def export_cooccurs_in_literature(self) -> None:
        # ---- Compound-Compound
        rows_cc = []
        for fp in self.cfg.linkdb_compound_compound_files:
            d = pd.read_csv(fp)
            for _, r in d.iterrows():
                id1 = _parse_py_dict(r.get("ID_1"))
                id2 = _parse_py_dict(r.get("ID_2"))
                c1 = _safe_int(id1.get("CID"))
                c2 = _safe_int(id2.get("CID"))
                if c1 is None or c2 is None or c1 == c2:
                    continue

                # de-dup as undirected (store ordered)
                a, b = (c1, c2) if c1 < c2 else (c2, c1)

                props = _parse_linkdb_block(r.get("Evidence"), key="ChemicalNeighbor")
                if not props:
                    continue

                rows_cc.append({
                    "CID": a,
                    "NeighborCID": b,
                    **props
                })

        if rows_cc:
            df_cc = pd.DataFrame(rows_cc)

            # aggregate duplicates across many CID_*.csv files
            agg = df_cc.groupby(["CID", "NeighborCID"], as_index=False).agg({
                "cooccurrence_score": "max",
                "article_count": "max",
                "pmids": _union_semicolon_lists,
                "dois": _union_semicolon_lists,
                "source": "first",
            })
            agg.to_csv(self.out_dir / "Compound_COOCCURS_IN_LITERATURE_Compound.csv", index=False)

        # ---- GeneSymbol-Compound
        rows_gc = []
        for fp in self.cfg.linkdb_gene_symbol_files:
            d = pd.read_csv(fp)
            for _, r in d.iterrows():
                id1 = _parse_py_dict(r.get("ID_1"))
                id2 = _parse_py_dict(r.get("ID_2"))
                sym = _safe_str(id1.get("GeneSymbol"))
                cid = _safe_int(id2.get("CID"))

                if not sym or sym.lower() == "nan" or cid is None:
                    continue

                props = _parse_linkdb_block(r.get("Evidence"), key="GeneSymbolChemicalNeighbor")
                if not props:
                    continue

                rows_gc.append({
                    "GeneSymbol": sym,
                    "CID": cid,
                    **props
                })

        if rows_gc:
            df_gc = pd.DataFrame(rows_gc)
            agg = df_gc.groupby(["GeneSymbol", "CID"], as_index=False).agg({
                "cooccurrence_score": "max",
                "article_count": "max",
                "pmids": _union_semicolon_lists,
                "dois": _union_semicolon_lists,
                "source": "first",
            })
            agg.to_csv(self.out_dir / "GeneSymbol_COOCCURS_IN_LITERATURE_Compound.csv", index=False)

    # -----------------------------
    # NEW: (:Compound)-[:TARGETS_GENE]->(:Gene)
    # -----------------------------
    def export_targets_gene(self) -> None:
        frames = []
        for fp in self.cfg.sdq_consolidatedcompoundtarget_files:
            d = pd.read_csv(fp)

            if "cid" not in d.columns or "geneid" not in d.columns:
                continue

            d = d.copy()
            d["CID"] = d["cid"].apply(_safe_int)
            d["GeneID"] = d["geneid"].apply(_safe_int)
            d = d.dropna(subset=["CID", "GeneID"])

            # keep SDQ fields as relationship properties
            d["source"] = "PubChemSDQ"
            d["collection"] = "consolidatedcompoundtarget"

            # drop raw duplicates to avoid confusion
            d = d.drop(columns=[c for c in ["cid", "geneid"] if c in d.columns])

            frames.append(d)

        if not frames:
            return

        df = pd.concat(frames, ignore_index=True)

        # Aggregate to 1 relationship per (CID, GeneID) (recommended)
        # - keep union of text fields, max of numeric-ish counts where meaningful
        text_cols = [c for c in df.columns if c not in {"CID", "GeneID"}]
        agg_map = {}
        for c in text_cols:
            if c in {"id"}:
                agg_map[c] = "max"
            else:
                agg_map[c] = _union_semicolon_lists

        out = df.groupby(["CID", "GeneID"], as_index=False).agg(agg_map)
        out.to_csv(self.out_dir / "Compound_TARGETS_GENE.csv", index=False)

    # -----------------------------
    # NEW: Metabolism normalization
    #   (:Compound)-[:SUBSTRATE_OF]->(:MetabolicTransformation)
    #   (:MetabolicTransformation)-[:PRODUCT]->(:Compound)
    #   (:MetabolicTransformation)-[:CATALYZED_BY]->(:Gene)
    #   (:MetabolicTransformation)-[:SUPPORTED_BY]->(:Publication)
    # -----------------------------
    def export_metabolism_normalized(self) -> None:
        mt = pd.read_csv(self.cfg.compound_transformation).copy()

        mt["SubstrateCID"] = mt["substratecid"].apply(_safe_int)
        mt["MetaboliteCID"] = mt["metabolitecid"].apply(_safe_int)
        mt["metconversion"] = mt["metconversion"].apply(_safe_str)
        mt["GeneID"] = mt["geneids"].apply(_safe_int)

        mt = mt.dropna(subset=["SubstrateCID", "MetaboliteCID"])

        # Stable event id
        mt["TransformationID"] = (
            mt["SubstrateCID"].astype(int).astype(str)
            + "->"
            + mt["MetaboliteCID"].astype(int).astype(str)
            + ":"
            + mt["metconversion"].fillna("")
        )

        # Nodes: MetabolicTransformation
        nodes = mt[["TransformationID", "metconversion"]].drop_duplicates().copy()
        nodes["source"] = "PubChemSDQ"
        nodes["collection"] = "chemblmetabolism"
        nodes.to_csv(self.out_dir / "MetabolicTransformation_NODES.csv", index=False)

        # Edge: Compound -SUBSTRATE_OF-> Transformation
        mt[["SubstrateCID", "TransformationID"]].drop_duplicates().rename(
            columns={"SubstrateCID": "CID"}
        ).to_csv(self.out_dir / "Compound_SUBSTRATE_OF_Transformation.csv", index=False)

        # Edge: Transformation -PRODUCT-> Compound
        mt[["TransformationID", "MetaboliteCID"]].drop_duplicates().rename(
            columns={"MetaboliteCID": "CID"}
        ).to_csv(self.out_dir / "Transformation_PRODUCT_Compound.csv", index=False)

        # Edge: Transformation -CATALYZED_BY-> Gene
        cat = mt.dropna(subset=["GeneID"])[["TransformationID", "GeneID"]].drop_duplicates()
        cat.to_csv(self.out_dir / "Transformation_CATALYZED_BY_Gene.csv", index=False)

        # Edge: Transformation -SUPPORTED_BY-> Publication (explode pmids/dois)
        pub_rows = []
        for _, r in mt[["TransformationID", "pmids", "dois"]].iterrows():
            tid = r["TransformationID"]

            for pmid in _split_pmids(r.get("pmids")):
                pub_rows.append((tid, f"PMID:{pmid}"))

            for doi in _split_dois(r.get("dois")):
                pub_rows.append((tid, f"DOI:{doi}"))

        if pub_rows:
            pubs = pd.DataFrame(pub_rows, columns=["TransformationID", "PublicationID"]).drop_duplicates()
            pubs.to_csv(self.out_dir / "Transformation_SUPPORTED_BY_Publication.csv", index=False)

    # -----------------------------
    # Run all exports
    # -----------------------------
    def export_all(self) -> None:
        self.export_core_from_all_data()
        self.export_compound_protein_edges()
        self.export_cooccurs_in_literature()
        self.export_targets_gene()
        self.export_metabolism_normalized()


# -----------------------------
# Example usage (call from your code)
# -----------------------------
def build_config_from_paths(
    all_data_connected: str,
    compound_gene_relationship: str,
    linkdb_compound_compound_glob: str,
    linkdb_gene_symbol_glob: str,
    sdq_glob: str,
    compound_transformation: str,
    out_dir: str,
) -> SchemaExportConfig:
    return SchemaExportConfig(
        all_data_connected=all_data_connected,
        compound_gene_relationship=compound_gene_relationship,
        linkdb_compound_compound_files=sorted(str(p) for p in Path().glob(linkdb_compound_compound_glob)),
        linkdb_gene_symbol_files=sorted(str(p) for p in Path().glob(linkdb_gene_symbol_glob)),
        sdq_consolidatedcompoundtarget_files=sorted(str(p) for p in Path().glob(sdq_glob)),
        compound_transformation=compound_transformation,
        out_dir=out_dir,
    )
