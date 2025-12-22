"""
NodesOntologyEnricher

This module defines a NodesOntologyEnricher class that:
    - reads the processed node CSVs produced by ChemGraphBuilder
    - calls external ontology/annotation services (MyGene.info, UniProt, PubChem, OLS)
    - adds ontology IDs / labels to each node table
    - writes enriched CSVs back to disk for Neo4j loading

Dependencies (install in your environment):

    pip install mygene bioservices pubchempy requests pandas

External services used:

    * MyGene.info (gene annotation)
      - https://mygene.info
      - Python wrapper: 'mygene'  :contentReference[oaicite:1]{index=1}

    * UniProt REST via bioservices.UniProt (protein + GO/InterPro)
      - https://www.uniprot.org  :contentReference[oaicite:2]{index=2}
      - https://bioservices.readthedocs.io/en/main/  :contentReference[oaicite:3]{index=3}

    * PubChemPy (to fetch cross-references such as ChEBI)
      - https://pubchempy.readthedocs.io  :contentReference[oaicite:4]{index=4}

    * EBI Ontology Lookup Service (OLS) for BAO/ChEBI term lookup
      - https://www.ebi.ac.uk/ols4/api  :contentReference[oaicite:5]{index=5}
"""

from __future__ import annotations

from collections import Counter
import datetime
import logging
import re
import time
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union
from biothings_client import get_client
import pandas as pd
import requests
from difflib import SequenceMatcher

try:
    import mygene
except ImportError:  # pragma: no cover
    mygene = None

try:
    from bioservices import UniProt
except ImportError:  # pragma: no cover
    UniProt = None

try:
    import pubchempy as pcp
except ImportError:  # pragma: no cover
    pcp = None


logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class NodesOntologyEnricher:
    """
    Enrich ChemGraphBuilder node tables with ontology identifiers
    using external web APIs / packages.

    Supported node types & ontologies (default):

        * Gene:
            - HGNC ID / symbol
            - Ensembl Gene ID
            - UniProt accession (if available via MyGene.info)
            - Selected GO terms (BP/MF/CC as joined lists)

        * Protein:
            - UniProt accession (from RefSeq via UniProt mapping)
            - GO term IDs (all)
            - InterPro accessions (IPR terms)

        * Compound:
            - ChEBI ID (CHEBI:xxxx)
            - MeSH ID(s)
            - ChEMBL ID(s)

        * BioAssay:
            - BAO term (best OLS match on assay name / activity name)
              (BAO:xxxx label & IRI)

    Assumed input files (your 'processed' node CSVs):

        data_dir/
            Gene_Properties_Processed.csv
            Protein_Properties_Processed.csv
            Compound_Properties_Processed.csv
            Assay_Properties_Processed.csv

    Output files (same directory, suffixed with '_WithOntologies'):

        Gene_Properties_WithOntologies.csv
        Protein_Properties_WithOntologies.csv
        Compound_Properties_WithOntologies.csv
        Assay_Properties_WithOntologies.csv
    """

    MYGENE_BATCH_SIZE = 1000
    UNIPROT_BATCH_SIZE = 500
    OLS_SLEEP_SECONDS = 0.2  # be polite to EBI OLS

    def __init__(
        self,
        data_dir: str = "Data/Nodes",
        species: str = "human",
        mygene_email: Optional[str] = None,
        ols_base_url: str = "https://www.ebi.ac.uk/ols4/api",
        http_timeout: int = 30,
    ):
        self._nlp = None   # lazy-initialised spaCy model
        self.data_dir = data_dir.rstrip("/")
        self.species = species
        self.ols_base_url = ols_base_url
        self.http_timeout = http_timeout

        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "ChemGraphBuilder-NodesOntologyEnricher/1.0"})

        # MyGene.info client
        if mygene is None:
            logger.warning(
                "mygene package not installed; gene ontology enrichment will be disabled."
            )
            self.mg = None
        else:
            self.mg = mygene.MyGeneInfo()
            if mygene_email:
                # not required but recommended by MyGene.info
                self.mg.email = mygene_email

        # UniProt client
        if UniProt is None:
            logger.warning(
                "bioservices.UniProt not installed; protein UniProt/GO/InterPro enrichment will be disabled."
            )
            self.uniprot = None
        else:
            self.uniprot = UniProt(verbose=False)

        if pcp is None:
            logger.warning(
                "pubchempy not installed; compound ChEBI/MeSH/ChEMBL enrichment will be disabled."
            )

    # ------------------------------------------------------------------
    # Public orchestration
    # ------------------------------------------------------------------
    def enrich_all(self) -> None:
        """Enrich all node types for which input files exist."""
        self.enrich_genes()
        self.enrich_proteins()
        self.enrich_compounds()
        self.enrich_assays()


    # ------------------------------------------------------------------
    # Compound enrichment via MyChem.info and PubChem
    # ------------------------------------------------------------------
    @staticmethod
    def _uniq(xs: List[str]) -> List[str]:
        seen = set()
        out = []
        for x in xs:
            if x and x not in seen:
                seen.add(x)
                out.append(x)
        return out

    @staticmethod
    def _count_non_empty_column(series: pd.Series) -> int:
        """
        Count rows where the ontology column is not empty.
        Treat NaN and "" as empty.
        """
        if series is None:
            return 0
        s = series.copy()
        mask = s.notna() & s.astype(str).str.strip().ne("")
        return int(mask.sum())


    def _save_run_summary(self, summary: dict, output_json_path: str):
        """
        Save ontology enrichment metadata for transparency.
        """
        import json
        with open(output_json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        logger.info("Saved ontology run summary to %s", output_json_path)

    def fetch_compound_xrefs_mychem_biothings(self, inchikeys: List[str], batch_size: int = 200
    ) -> Dict[str, Dict[str, List[str]]]:
        """
        Retrieve all ontology/xref data using MyChem (BioThings client),
        including PubChem cross-references.
        """

        chem = get_client("chem")

        keys = sorted({ik.strip() for ik in inchikeys if isinstance(ik, str) and ik.strip()})
        results: Dict[str, Dict[str, List[str]]] = {}

        if not keys:
            return results

        fields = ",".join([
            "chebi",
            "chembl",
            "drugbank",
            "unii",
            "pubchem.xrefs",        # <-- unlocks: pubmed, gene, mim, taxonomy, mesh, kegg, hmdb, dtxsid, cas, patent, url
        ])

        for i in range(0, len(keys), batch_size):
            batch = keys[i:i+batch_size]

            logger.info(f"Querying MyChem batch {i+1}-{i+len(batch)} …")

            try:
                docs = chem.getchems(batch, fields=fields, as_dataframe=False)
            except Exception as e:
                logger.error("MyChem batch failed: %s", e)
                continue

            for doc in docs:
                if not isinstance(doc, dict):
                    continue

                ik = doc.get("_id") or doc.get("query")
                if not ik:
                    continue

                out = {}

                # -------------------------------
                # Direct ontologies (ChEBI/ChEMBL/DrugBank/UNII)
                # -------------------------------
                def collect(source, keys):
                    vals = []
                    if isinstance(source, dict):
                        for k in keys:
                            v = source.get(k)
                            if isinstance(v, str):
                                vals.append(v)
                            elif isinstance(v, list):
                                vals.extend(v)
                    elif isinstance(source, list):
                        for item in source:
                            if isinstance(item, dict):
                                for k in keys:
                                    v = item.get(k)
                                    if isinstance(v, str):
                                        vals.append(v)
                                    elif isinstance(v, list):
                                        vals.extend(v)
                    return vals

                # ChEBI
                out["ChEBI"] = self._uniq(collect(doc.get("chebi"), ["chebi_id", "id"]))

                # ChEMBL
                out["ChEMBL"] = self._uniq(collect(doc.get("chembl"),
                                            ["molecule_chembl_id", "drug_chembl_id", "id"]))

                # DrugBank
                out["DrugBank"] = self._uniq(collect(doc.get("drugbank"), ["id"]))

                # UNII + CAS
                unii_block = doc.get("unii")
                out["UNII"] = self._uniq(collect(unii_block, ["unii", "id"]))
                out["CAS_RN"] = self._uniq(collect(unii_block, ["registry_number"]))

                # -------------------------------
                # PubChem xrefs: full ontology set
                # -------------------------------
                xref = doc.get("pubchem", {}).get("xrefs", {})

                def gx(field):
                    v = xref.get(field, [])
                    if isinstance(v, str):
                        return [v]
                    if isinstance(v, list):
                        return [str(x) for x in v]
                    return []

                out["PubMedID"]   = self._uniq(gx("pubmed"))
                out["PatentID"]   = self._uniq(gx("patent"))
                out["MIMID"]      = self._uniq(gx("mim"))
                out["GeneID"]     = self._uniq(gx("gene"))
                out["TaxonomyID"] = self._uniq(gx("taxonomy"))
                # out["MeSH"]       = self._uniq(gx("mesh"))
                # out["KEGG"]       = self._uniq(gx("kegg"))
                # out["HMDB"]       = self._uniq(gx("hmdb"))
                # out["DTXSID"]     = self._uniq(gx("dtxsid"))
                out["DBURL"]      = self._uniq(gx("url"))
                out["CAS_RN"]    += self._uniq(gx("cas"))  # add CAS from pubchem.xrefs

                results[ik] = out

        return results

    def _fetch_compound_xrefs_pubchem_batch(
        self,
        cids: List[int],
        batch_size: int = 100,
        max_retries: int = 3,
        retry_sleep: float = 5.0,
    ) -> Dict[int, Dict[str, List[str]]]:
        """
        Fetch PubChem xrefs in batch mode by CID using PUG REST.

        Processes up to `batch_size` CIDs per HTTP call and retries on failures.

        Returns:
            {
            cid: {
                "RegistryID": [...],
                "PubMedID": [...],
                "MIMID": [...],
                "GeneID": [...],
                "PatentID": [...],
                "DBURL": [...],
                "TaxonomyID": [...],
            },
            ...
            }
        """
        results: Dict[int, Dict[str, List[str]]] = {}
        uniq_cids = sorted({int(c) for c in cids if pd.notna(c)})

        if not uniq_cids:
            return results

        xref_types = "RegistryID,PubMedID,MIMID,GeneID,PatentID,DBURL,TaxonomyID"
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"

        logger.info(
            "Fetching PubChem xrefs for %d unique CIDs (batch_size=%d)",
            len(uniq_cids),
            batch_size,
        )

        for i in range(0, len(uniq_cids), batch_size):
            batch = uniq_cids[i : i + batch_size]
            payload = {"cid": ",".join(str(c) for c in batch)}

            attempt = 0
            success = False

            while attempt < max_retries and not success:
                attempt += 1
                logger.info(
                    "PubChem xrefs batch %d-%d (size=%d), attempt %d/%d",
                    i + 1,
                    i + len(batch),
                    len(batch),
                    attempt,
                    max_retries,
                )

                try:
                    resp = self.session.post(
                        f"{base_url}/xrefs/{xref_types}/JSON",
                        data=payload,
                        timeout=self.http_timeout,
                    )

                    # 404: no xrefs – treat as empty batch and stop retrying
                    if resp.status_code == 404:
                        logger.warning(
                            "PubChem xrefs batch %d-%d returned 404 (no data).",
                            i + 1,
                            i + len(batch),
                        )
                        success = True
                        data = {}
                        break

                    # 504 or other server error – retry
                    if resp.status_code == 504:
                        logger.warning(
                            "PubChem xrefs batch %d-%d got 504 (timeout).",
                            i + 1,
                            i + len(batch),
                        )
                        if attempt < max_retries:
                            time.sleep(retry_sleep)
                        continue

                    resp.raise_for_status()
                    data = resp.json()
                    success = True

                except Exception as e:
                    logger.error(
                        "PubChem xrefs batch %d-%d failed on attempt %d/%d: %s",
                        i + 1,
                        i + len(batch),
                        attempt,
                        max_retries,
                        e,
                    )
                    if attempt < max_retries:
                        time.sleep(retry_sleep)

            if not success:
                logger.error(
                    "PubChem xrefs batch %d-%d permanently failed after %d attempts; "
                    "these CIDs will have no PubChem xrefs.",
                    i + 1,
                    i + len(batch),
                    max_retries,
                )
                continue

            info_list = (
                data.get("InformationList", {})
                .get("Information", [])
                if isinstance(data, dict)
                else []
            )

            for info in info_list:
                if not isinstance(info, dict):
                    continue

                cid_val = info.get("CID")
                if cid_val is None:
                    continue

                cid_val = int(cid_val)
                info = {k: v for k, v in info.items() if k != "CID"}

                out: Dict[str, List[str]] = {}
                for k, v in info.items():
                    if v is None:
                        continue
                    if isinstance(v, list):
                        vals = [str(x) for x in v if x is not None]
                    else:
                        vals = [str(v)]
                    if vals:
                        out[k] = vals

                if out:
                    results[cid_val] = out
            time.sleep(0.5)  # be polite to PubChem
        return results

    def enrich_compounds(
        self,
        input_name: str = "Compound_Properties_Processed.csv",
        output_name: str = "Compound_Properties_WithOntologies.csv",
        inchikey_col: str = "InChIKey",
        cid_col: str = "CompoundID",
    ) -> Optional[pd.DataFrame]:

        input_path = f"{self.data_dir}/{input_name}"

        try:
            df = pd.read_csv(input_path, low_memory=False)
        except FileNotFoundError:
            logger.error("Compound file not found: %s", input_path)
            return None

        logger.info("Running compound ontology enrichment (MyChem + PubChem batch)…")

        # ------------------------------------------------------
        # 1) Fetch all ontologies via MyChem BioThings
        # ------------------------------------------------------
        if inchikey_col not in df.columns:
            logger.error("No InChIKey column '%s' found", inchikey_col)
            return None

        inchikeys = (
            df[inchikey_col].dropna().astype(str).str.strip().tolist()
        )
        uniq_iks = sorted(set(inchikeys))

        logger.info("Enriching %d unique compounds via MyChem…", len(uniq_iks))
        mychem_xrefs = self.fetch_compound_xrefs_mychem_biothings(uniq_iks)

        # ------------------------------------------------------
        # 2) Fetch xrefs via PubChem PUG REST (batch by CID)
        # ------------------------------------------------------
        pubchem_xrefs: Dict[int, Dict[str, List[str]]] = {}
        cid_series = None

        if cid_col in df.columns:
            cid_series = df[cid_col]
        elif "CID" in df.columns:
            logger.warning(
                "CID column '%s' not found; falling back to 'CID'.",
                cid_col,
            )
            cid_col = "CID"
            cid_series = df[cid_col]
        else:
            logger.warning(
                "No CID column ('%s' or 'CID') found; PubChem xrefs will be skipped.",
                cid_col,
            )

        if cid_series is not None:
            cid_values = cid_series.dropna().astype(int).tolist()
            pubchem_xrefs = self._fetch_compound_xrefs_pubchem_batch(
                cid_values,
                batch_size=100,
            )

        # ------------------------------------------------------
        # 3) Build output columns (IDs + per-ontology sources)
        # ------------------------------------------------------
        # logical_key = internal name we use across MyChem/PubChem
        output_fields: Dict[str, str] = {
            "ChEBI_IDs": "ChEBI",
            "ChEMBL_IDs": "ChEMBL",
            "DrugBank_IDs": "DrugBank",
            "UNII_IDs": "UNII",
            "CAS_RN": "CAS_RN",
            "OMIM_IDs": "MIMID",
            "GeneIDs": "GeneID",
            "PubMed_IDs": "PubMedID",
            "Patent_IDs": "PatentID",
            "Taxonomy_IDs": "TaxonomyID",
            "DBURLs": "DBURL",
        }

        # Map ontology ID columns → source columns
        source_columns: Dict[str, str] = {
            "ChEBI_IDs": "ChEBI_Sources",
            "ChEMBL_IDs": "ChEMBL_Sources",
            "DrugBank_IDs": "DrugBank_Sources",
            "UNII_IDs": "UNII_Sources",
            "CAS_RN": "CAS_RN_Sources",
            "OMIM_IDs": "OMIM_Sources",
            "GeneIDs": "GeneIDs_Sources",
            "PubMed_IDs": "PubMed_Sources",
            "Patent_IDs": "Patent_Sources",
            "Taxonomy_IDs": "Taxonomy_Sources",
            "DBURLs": "DBURLs_Sources",
        }

        # Track per-ontology, per-source counts for the run summary
        ontology_source_counts: Dict[str, Dict[str, int]] = {
            out_col: {"MyChem": 0, "PubChem": 0}
            for out_col in output_fields.keys()
        }

        def get_joined_with_sources(
            ik: str,
            cid: Optional[int],
            logical_key: str,
        ) -> Tuple[str, str]:
            """
            Merge MyChem + PubChem values for a logical ontology key
            and return (joined_ids, sources_string).

            sources_string is in {"", "MyChem", "PubChem", "MyChem;PubChem"}.
            """
            vals: List[str] = []
            sources: Set[str] = set()

            # MyChem values
            if ik:
                m_vals = mychem_xrefs.get(ik, {}).get(logical_key, [])
                if m_vals:
                    vals.extend(m_vals)
                    sources.add("MyChem")

            # PubChem values – only for selected keys
            if cid is not None and cid in pubchem_xrefs:
                px = pubchem_xrefs[cid]

                if logical_key in {"MIMID", "GeneID", "PubMedID", "PatentID", "TaxonomyID", "DBURL"}:
                    p_vals = px.get(logical_key, [])
                    if p_vals:
                        vals.extend(p_vals)
                        sources.add("PubChem")

                elif logical_key == "CAS_RN":
                    # Derive CAS from RegistryID if it matches digits-digits-digit pattern
                    regs = px.get("RegistryID", [])
                    cas_candidates: List[str] = []
                    for reg in regs:
                        if reg and "-" in reg and reg.replace("-", "").isdigit():
                            cas_candidates.append(reg)
                    if cas_candidates:
                        vals.extend(cas_candidates)
                        sources.add("PubChem")

            vals = [v for v in vals if v]
            joined_vals = "|".join(sorted(set(vals))) if vals else ""
            joined_sources = ";".join(sorted(sources)) if sources else ""

            return joined_vals, joined_sources

        # Prepare buffers for all ID + source columns
        id_columns_data: Dict[str, List[str]] = {out_col: [] for out_col in output_fields.keys()}
        source_columns_data: Dict[str, List[str]] = {src_col: [] for src_col in source_columns.values()}

        # Iterate once over rows and fill all ontology columns
        for _, row in df.iterrows():
            ik_val = str(row[inchikey_col]).strip() if pd.notna(row[inchikey_col]) else ""
            cid_val: Optional[int] = None
            if cid_series is not None:
                cid_raw = row.get(cid_col)
                if pd.notna(cid_raw):
                    try:
                        cid_val = int(cid_raw)
                    except Exception:
                        cid_val = None

            for out_col, logical_key in output_fields.items():
                ids_str, src_str = get_joined_with_sources(ik_val, cid_val, logical_key)

                id_columns_data[out_col].append(ids_str)

                src_col_name = source_columns[out_col]
                source_columns_data[src_col_name].append(src_str)

                # Update per-source counts
                if "MyChem" in src_str:
                    ontology_source_counts[out_col]["MyChem"] += 1
                if "PubChem" in src_str:
                    ontology_source_counts[out_col]["PubChem"] += 1

        # Attach ID + source columns to the DataFrame
        for out_col, values in id_columns_data.items():
            df[out_col] = values
        for src_col, values in source_columns_data.items():
            df[src_col] = values

        # ------------------------------------------------------
        # 4) Ontology counts summary + transparency metadata
        # ------------------------------------------------------
        ontology_counts: Dict[str, int] = {}
        for out_col in output_fields.keys():
            count = self._count_non_empty_column(df[out_col])
            ontology_counts[out_col] = int(count)

        num_total = len(uniq_iks)

        # Percentage coverage for each ontology (any source)
        ontology_percent = {
            k: (v / num_total * 100.0 if num_total > 0 else 0.0)
            for k, v in ontology_counts.items()
        }

        # Compounds with no ontology fields at all (based on final DF)
        missing_all: List[str] = []
        for _, row in df.iterrows():
            ik_val = str(row[inchikey_col]).strip() if pd.notna(row[inchikey_col]) else ""
            has_any = any(
                bool(str(row[col]).strip())
                for col in output_fields.keys()
            )
            if not has_any and ik_val:
                missing_all.append(ik_val)

        num_missing = len(missing_all)

        logger.info("=== Ontology Retrieval Summary ===")
        logger.info("Total compounds processed: %d", num_total)
        logger.info("Successful hits (any ontology): %d", num_total - num_missing)
        logger.info("Missing entries (no ontologies at all): %d", num_missing)
        logger.info("Ontology counts / coverage (all sources):")
        for k, v in ontology_counts.items():
            logger.info(
                "  %s: %d  (%.2f%%) [MyChem=%d, PubChem=%d]",
                k,
                v,
                ontology_percent[k],
                ontology_source_counts[k]["MyChem"],
                ontology_source_counts[k]["PubChem"],
            )

        # Build final summary dict
        run_summary = {
            "total_compounds_processed": int(num_total),
            "total_unique_inchikeys": int(num_total),
            "total_successful_hits": int(num_total - num_missing),
            "total_missing_entries": int(num_missing),
            "missing_inchikeys": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "ontology_source_counts": ontology_source_counts,  # <-- NEW
            "retrieved_fields": list(output_fields.keys()),
            "input_file": input_path,
            "output_file": f"{self.data_dir}/{output_name}",
            "inchikey_column": inchikey_col,
            "cid_column": cid_col if cid_series is not None else None,
            "mychem_fields": ["chebi", "chembl", "drugbank", "unii", "pubchem.xrefs"],
            "pubchem_xref_types": [
                "RegistryID",
                "PubMedID",
                "GeneID",
                "PatentID",
                "DBURL",
                "TaxonomyID",
            ],
            "created_at": datetime.datetime.utcnow().isoformat() + "Z",
            "batch_size_mychem": 200,
            "batch_size_pubchem": 100,
            "source": "MyChem.info(BioThings);PubChem PUG REST",
        }

        summary_path = f"{self.data_dir}/Compound_Ontology_RunSummary.json"
        self._save_run_summary(run_summary, summary_path)

        # ------------------------------------------------------
        # 5) Save enriched version
        # ------------------------------------------------------
        df["OntologyEnrichedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.MyChem+PubChem/1.1"
        df["OntologySources"] = "MyChem.info;PubChem"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved enriched compound table to %s", out_path)

        return df


    # ------------------------------------------------------------------
    # Gene enrichment via MyGene.info (BioThings client)
    # ------------------------------------------------------------------
    def fetch_gene_xrefs_mygene(
        self,
        gene_ids: List[Union[int, str]],
        batch_size: int = 1000,
        species: str = "human",
    ) -> Dict[str, Dict[str, List[str]]]:
        """
        Retrieve gene-level ontologies / xrefs using MyGene (BioThings client).

        Returned structure:
            {
              "7157": {
                  "HGNC_Symbol": [...],
                  "HGNC_ID": [...],
                  "EnsemblGene": [...],
                  "GO_BP": [...],      # GO:0008150 etc
                  "GO_MF": [...],
                  "GO_CC": [...],
                  "Reactome": [...],   # R-HSA-xxxx
                  "OMIM": [...],
                  # DOID/MONDO/HPO are *not* currently populated here
              },
              ...
            }
        """

        mg = get_client("gene")

        # Normalise IDs to strings
        ids = sorted(
            {
                str(g).strip()
                for g in gene_ids
                if g is not None and not (isinstance(g, float) and pd.isna(g))
            }
        )
        results: Dict[str, Dict[str, List[str]]] = {}

        if not ids:
            return results

        # Fields we actually care about
        # (DOID/MONDO/HPO are not reliably exposed per gene via MyGene alone.)
        fields = ",".join(
            [
                "symbol",
                "HGNC",           # often present as "HGNC"
                "hgnc",           # sometimes lowercase
                "ensembl",
                "go",
                "pathway.reactome",
                "omim",
            ]
        )

        logger.info(
            "Fetching gene ontologies from MyGene.info for %d unique GeneIDs "
            "(batch_size=%d, species=%s)",
            len(ids),
            batch_size,
            species,
        )

        for i in range(0, len(ids), batch_size):
            batch = ids[i : i + batch_size]
            logger.info(
                "MyGene batch %d-%d (size=%d)",
                i + 1,
                i + len(batch),
                len(batch),
            )

            try:
                docs = mg.getgenes(
                    batch,
                    fields=fields,
                    species=species,
                    as_dataframe=False,
                )
            except Exception as e:
                logger.error("MyGene batch %d-%d failed: %s", i + 1, i + len(batch), e)
                continue

            for doc in docs:
                if not isinstance(doc, dict):
                    continue

                gid = doc.get("_id") or doc.get("entrezgene") or doc.get("query")
                if gid is None:
                    continue
                gid_str = str(gid)

                out: Dict[str, List[str]] = {}

                # -------------------------------
                # HGNC Symbol
                # -------------------------------
                symbol = doc.get("symbol")
                if isinstance(symbol, str) and symbol.strip():
                    out["HGNC_Symbol"] = [symbol.strip()]

                # -------------------------------
                # HGNC ID
                # -------------------------------
                hgnc_ids: List[str] = []
                for key in ("HGNC", "hgnc", "HGNC_ID"):
                    v = doc.get(key)
                    if isinstance(v, str):
                        hgnc_ids.append(v)
                    elif isinstance(v, (int, float)):
                        hgnc_ids.append(str(v))
                    elif isinstance(v, list):
                        hgnc_ids.extend(str(x) for x in v if x)
                if hgnc_ids:
                    out["HGNC_ID"] = self._uniq(hgnc_ids)

                # -------------------------------
                # Ensembl gene IDs
                # -------------------------------
                ens_ids: List[str] = []
                ens = doc.get("ensembl")
                if isinstance(ens, dict):
                    v = (
                        ens.get("gene")
                        or ens.get("geneid")
                        or ens.get("id")
                    )
                    if isinstance(v, str):
                        ens_ids.append(v)
                    elif isinstance(v, list):
                        ens_ids.extend(str(x) for x in v if x)
                elif isinstance(ens, list):
                    for item in ens:
                        if not isinstance(item, dict):
                            continue
                        v = (
                            item.get("gene")
                            or item.get("geneid")
                            or item.get("id")
                        )
                        if isinstance(v, str):
                            ens_ids.append(v)
                        elif isinstance(v, list):
                            ens_ids.extend(str(x) for x in v if x)
                if ens_ids:
                    out["EnsemblGene"] = self._uniq(ens_ids)

                # -------------------------------
                # GO: BP / MF / CC (GO IDs)
                # -------------------------------
                go = doc.get("go", {})

                def collect_go(cat: str) -> List[str]:
                    vals: List[str] = []
                    sub = go.get(cat, [])
                    if isinstance(sub, dict):
                        sub = [sub]
                    if isinstance(sub, list):
                        for term in sub:
                            if not isinstance(term, dict):
                                continue
                            tid = term.get("id")
                            if tid:
                                vals.append(str(tid))
                    return self._uniq(vals) if vals else []

                bp_ids = collect_go("BP")
                mf_ids = collect_go("MF")
                cc_ids = collect_go("CC")

                if bp_ids:
                    out["GO_BP"] = bp_ids
                if mf_ids:
                    out["GO_MF"] = mf_ids
                if cc_ids:
                    out["GO_CC"] = cc_ids

                # -------------------------------
                # Reactome pathways (IDs)
                # -------------------------------
                react_ids: List[str] = []
                pathway = doc.get("pathway") or {}
                # MyGene uses pathway.reactome
                reactome = None
                if isinstance(pathway, dict):
                    reactome = pathway.get("reactome")

                if isinstance(reactome, dict):
                    reactome = [reactome]

                if isinstance(reactome, list):
                    for pw in reactome:
                        if not isinstance(pw, dict):
                            continue
                        pid = pw.get("id")
                        if pid:
                            react_ids.append(str(pid))

                if react_ids:
                    out["Reactome"] = self._uniq(react_ids)

                # -------------------------------
                # OMIM IDs
                # -------------------------------
                omim_ids: List[str] = []
                omim = doc.get("omim")
                if isinstance(omim, (str, int)):
                    omim_ids.append(str(omim))
                elif isinstance(omim, list):
                    for o in omim:
                        if isinstance(o, (str, int)):
                            omim_ids.append(str(o))
                        elif isinstance(o, dict):
                            oid = o.get("id")
                            if isinstance(oid, (str, int)):
                                omim_ids.append(str(oid))
                if omim_ids:
                    out["OMIM"] = self._uniq(omim_ids)

                # NOTE: DOID / MONDO / HPO are *not* collected here.
                # They typically require disease-centric KPs (Monarch, MyDisease, etc.).

                if out:
                    results[gid_str] = out

        return results


    # ------------------------------------------------------------------
    # Gene enrichment via MyGene.info
    # ------------------------------------------------------------------
    def _batch_iter(self, items: List, batch_size: int) -> Iterable[List]:
        for i in range(0, len(items), batch_size):
            yield items[i : i + batch_size]

    def enrich_genes(
        self,
        input_name: str = "Gene_Properties_Processed.csv",
        output_name: str = "Gene_Properties_WithOntologies.csv",
        gene_id_col: str = "GeneID",
    ) -> Optional[pd.DataFrame]:
        """
        Enrich Gene nodes with ontology cross-references using MyGene.info.

        Target ontology columns in the output CSV:

            GeneID                     (input)
            HGNC_Symbol                (MyGene)
            HGNC_IDs                   (MyGene)
            EnsemblGene_IDs            (MyGene)
            GO_BP_IDs                  (MyGene: go.BP IDs)
            GO_MF_IDs                  (MyGene: go.MF IDs)
            GO_CC_IDs                  (MyGene: go.CC IDs)
            Reactome_Pathway_IDs       (MyGene: pathway.reactome IDs)
            OMIM_IDs                   (MyGene: omim IDs)

            DOID_IDs                   (reserved, currently not populated)
            MONDO_IDs                  (reserved, currently not populated)
            HPO_Terms                  (reserved, currently not populated)

        A detailed JSON run summary is written to:

            Data/Nodes/Gene_Ontology_RunSummary.json
        """

        input_path = f"{self.data_dir}/{input_name}"

        try:
            df = pd.read_csv(input_path, low_memory=False)
        except FileNotFoundError:
            logger.error("Gene file not found: %s", input_path)
            return None

        logger.info("Running gene ontology enrichment via MyGene.info…")

        # ------------------------------------------------------
        # 1) Validate gene ID column
        # ------------------------------------------------------
        if gene_id_col not in df.columns:
            logger.error("No GeneID column '%s' found in %s", gene_id_col, input_path)
            return None

        gene_ids = df[gene_id_col].dropna().tolist()
        uniq_gene_ids = sorted({str(g).strip() for g in gene_ids})

        logger.info(
            "Enriching %d unique genes via MyGene.info (by %s)…",
            len(uniq_gene_ids),
            gene_id_col,
        )

        gene_xrefs = self.fetch_gene_xrefs_mygene(uniq_gene_ids)

        # ------------------------------------------------------
        # 2) Build output columns
        # ------------------------------------------------------
        def get_joined(gid: str, field: str) -> str:
            vals = gene_xrefs.get(gid, {}).get(field, [])
            vals = [v for v in vals if v]
            return "|".join(sorted(set(vals))) if vals else ""

        # Logical keys in gene_xrefs -> output columns
        output_fields = {
            "HGNC_Symbol": "HGNC_Symbol",
            "HGNC_IDs": "HGNC_ID",
            "EnsemblGene_IDs": "EnsemblGene",
            "GO_BP_IDs": "GO_BP",
            "GO_MF_IDs": "GO_MF",
            "GO_CC_IDs": "GO_CC",
            "Reactome_Pathway_IDs": "Reactome",
            "OMIM_IDs": "OMIM",
            # Reserve disease/phenotype ontologies
            # "DOID_IDs": "DOID",    # not populated (requires external KPs)
            # "MONDO_IDs": "MONDO",  # not populated
            # "HPO_Terms": "HPO",    # not populated
        }

        # Fill columns
        for out_col, logical_key in output_fields.items():
            values: List[str] = []
            for _, row in df.iterrows():
                gid_val = row.get(gene_id_col)
                gid_str = str(gid_val).strip() if pd.notna(gid_val) else ""
                if logical_key in {"DOID", "MONDO", "HPO"}:
                    # currently not populated from MyGene
                    values.append("")
                else:
                    values.append(get_joined(gid_str, logical_key))
            df[out_col] = values

        # ------------------------------------------------------
        # 3) Ontology counts summary + transparency metadata
        # ------------------------------------------------------
        ontology_counts: Dict[str, int] = {}
        for out_col in output_fields.keys():
            count = self._count_non_empty_column(df[out_col])
            ontology_counts[out_col] = int(count)

        num_total = len(uniq_gene_ids)

        ontology_percent = {
            k: (v / num_total * 100.0 if num_total > 0 else 0.0)
            for k, v in ontology_counts.items()
        }

        # Identify genes with no ontologies at all (based on final DF)
        missing_all: List[str] = []
        for _, row in df.iterrows():
            gid_val = row.get(gene_id_col)
            gid_str = str(gid_val).strip() if pd.notna(gid_val) else ""
            has_any = any(
                bool(str(row[col]).strip())
                for col in [
                    "HGNC_Symbol",
                    "HGNC_IDs",
                    "EnsemblGene_IDs",
                    "GO_BP_IDs",
                    "GO_MF_IDs",
                    "GO_CC_IDs",
                    "Reactome_Pathway_IDs",
                    "OMIM_IDs",
                ]
            )
            if not has_any and gid_str:
                missing_all.append(gid_str)

        num_missing = len(missing_all)

        logger.info("=== Gene Ontology Retrieval Summary ===")
        logger.info("Total genes processed: %d", num_total)
        logger.info(
            "Successful hits (any ontology): %d",
            num_total - num_missing,
        )
        logger.info(
            "Missing entries (no ontologies at all): %d",
            num_missing,
        )
        logger.info("Ontology counts / coverage:")
        for k, v in ontology_counts.items():
            logger.info("  %s: %d  (%.2f%%)", k, v, ontology_percent[k])

        # Per-ontology source mapping (for transparency)
        ontology_sources = {
            "HGNC_Symbol": "MyGene.info (symbol)",
            "HGNC_IDs": "MyGene.info (HGNC / hgnc)",
            "EnsemblGene_IDs": "MyGene.info (ensembl)",
            "GO_BP_IDs": "MyGene.info (go.BP)",
            "GO_MF_IDs": "MyGene.info (go.MF)",
            "GO_CC_IDs": "MyGene.info (go.CC)",
            "Reactome_Pathway_IDs": "MyGene.info (pathway.reactome)",
            "OMIM_IDs": "MyGene.info (omim)",
            # "DOID_IDs": "Not populated – requires disease-centric KPs (e.g. Monarch/MyDisease)",
            # "MONDO_IDs": "Not populated – requires disease-centric KPs (e.g. ClinGen/Monarch)",
            # "HPO_Terms": "Not populated – requires phenotype-centric KPs (e.g. HPO/Monarch)",
        }

        # Build and save JSON summary
        run_summary = {
            "total_genes_processed": int(num_total),
            "total_unique_gene_ids": int(num_total),
            "total_successful_hits": int(num_total - num_missing),
            "total_missing_entries": int(num_missing),
            "missing_gene_ids": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "retrieved_fields": list(output_fields.keys()),
            "ontology_sources": ontology_sources,
            "input_file": input_path,
            "output_file": f"{self.data_dir}/{output_name}",
            "gene_id_column": gene_id_col,
            "mygene_fields": [
                "symbol",
                "HGNC/hgnc",
                "ensembl",
                "go (BP/MF/CC)",
                "pathway.reactome",
                "omim",
            ],
            "created_at": datetime.datetime.utcnow().isoformat() + "Z",
            "batch_size_mygene": 1000,
            "source": "MyGene.info (via BioThings client)",
        }

        summary_path = f"{self.data_dir}/Gene_Ontology_RunSummary.json"
        self._save_run_summary(run_summary, summary_path)

        # ------------------------------------------------------
        # 4) Save enriched version
        # ------------------------------------------------------
        df["OntologyEnrichedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.GenesMyGene/1.0"
        df["OntologySources"] = "MyGene.info"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved enriched gene table to %s", out_path)

        return df


    # ------------------------------------------------------------------
    # Protein enrichment via UniProt (RefSeq -> UniProt -> GO / InterPro / Pfam / PDB / Reactome)
    # ------------------------------------------------------------------
    @staticmethod
    def _looks_like_uniprot_accession(acc: str) -> bool:
        """
        Heuristic check if a string looks like a UniProtKB accession.

        Supports 6-char accessions optionally followed by an isoform suffix,
        e.g. 'P08684' or 'P08684-2'.
        """
        import re

        if not acc:
            return False

        acc = acc.strip()
        core = acc.split("-", 1)[0]

        # Robust pattern: 6 alphanumeric chars (standard UniProt accessions).
        return bool(re.match(r"^[A-Z0-9]{6}$", core))

    def _fetch_uniprot_protein_annotations(
        self,
        accessions: List[str],
        batch_size: int = 200,
    ) -> Dict[str, Dict[str, List[str]]]:
        """
        Fetch GO (BP/MF/CC) and a broad set of ontology / pathway / family /
        disease / interaction cross-references for UniProt accessions
        via the UniProt REST 'stream' endpoint.

        Returns (per UniProt accession) a dict like:
            {
                "GO_BP": [...],
                "GO_MF": [...],
                "GO_CC": [...],
                "InterPro": [...],
                "Pfam": [...],
                "PDB": [...],
                "Reactome": [...],
                "KEGG": [...],
                "UniPathway": [...],
                "PlantReactome": [...],
                "PANTHER": [...],
                "PIRSF": [...],
                "PRINTS": [...],
                "PROSITE": [...],
                "Gene3D": [...],
                "FunFam": [...],
                "SFLD": [...],
                "SUPFAM": [...],
                "SMART": [...],
                "CDD": [...],
                "HAMAP": [...],
                "PRO": [...],
                "STRING": [...],
                "BioGRID": [...],
                "IntAct": [...],
                "NCBIGene": [...],
                "Ensembl": [...],
                "Orphanet": [...],
                "MIM": [...],
            }
        """
        annotations: Dict[str, Dict[str, List[str]]] = {}

        if not accessions:
            return annotations
        if self.session is None:
            logger.warning("No HTTP session available; skipping UniProt protein annotations.")
            return annotations

        # ------------------------------------------------------------------
        # Define which UniProt xref_* fields to treat as ontologies / xrefs
        # ------------------------------------------------------------------
        # Keys = UniProt column names (xref_*)
        # Values = labels used internally + for output mapping.
        xref_ontology_fields: Dict[str, str] = {
            # Domain / family / signature resources
            "xref_interpro": "InterPro",
            "xref_pfam": "Pfam",
            "xref_prosite": "PROSITE",
            "xref_panther": "PANTHER",
            "xref_pirsf": "PIRSF",
            "xref_prints": "PRINTS",
            "xref_gene3d": "Gene3D",
            "xref_funfam": "FunFam",         # CATH FunFams
            "xref_sfld": "SFLD",
            "xref_supfam": "SUPFAM",
            "xref_smart": "SMART",
            "xref_cdd": "CDD",
            "xref_hamap": "HAMAP",
            "xref_pro": "PRO",

            # Structure / pathways
            "xref_pdb": "PDB",
            "xref_reactome": "Reactome",
            "xref_kegg": "KEGG",
            "xref_unipathway": "UniPathway",
            "xref_plantreactome": "PlantReactome",

            # Interactions
            "xref_string": "STRING",
            "xref_biogrid": "BioGRID",
            "xref_intact": "IntAct",

            # Gene / disease identifiers
            "xref_geneid": "NCBIGene",
            "xref_ensembl": "Ensembl",
            "xref_mim": "MIM",              # OMIM
            "xref_orphanet": "Orphanet",
        }
        # If you want to add more later, just extend this dict using
        # the xref_* field names from UniProt's `return_fields` help. :contentReference[oaicite:2]{index=2}

        url = "https://rest.uniprot.org/uniprotkb/stream"

        # Order is important – we will parse by column position, not header text.
        # 0: accession
        # 1: go_p
        # 2: go_f
        # 3: go_c
        # 4..N: all selected xref_* fields in the order defined above.
        xref_field_list = list(xref_ontology_fields.keys())
        fields = "accession,go_p,go_f,go_c," + ",".join(xref_field_list)

        def _split_go_field(raw: str) -> List[str]:
            """
            Extract GO IDs (GO:########) from UniProt GO columns.

            Handles formats like:
              'amide metabolic process [GO:0043603]; cholesterol metabolic process [GO:0008203]; ...'
            possibly spanning multiple lines.

            We ignore labels/evidence and just keep the GO IDs.
            """
            if not raw:
                return []

            import re

            text = str(raw)
            ids = re.findall(r"GO:\d{7}", text)

            seen = set()
            out: List[str] = []
            for go_id in ids:
                if go_id not in seen:
                    seen.add(go_id)
                    out.append(go_id)
            return out

        def _split_db_field(raw: str) -> List[str]:
            """
            Generic parser for xref_* columns.

            UniProt xref_* fields are typically ';'-separated, with optional
            labels before ':' and descriptions after a space, e.g.:

                'InterPro:IPR001234 description; InterPro:IPR005678 other desc'

            We keep only the identifier part (after ':' if present), and
            take the first whitespace-separated token.
            """
            if not raw:
                return []
            entries: List[str] = []
            for piece in str(raw).split(";"):
                piece = piece.strip()
                if not piece:
                    continue
                # Split off label like 'InterPro:' or 'Pfam:'
                if ":" in piece:
                    piece = piece.split(":", 1)[1]
                # Keep first token (drop descriptions)
                piece = piece.split()[0]
                if piece:
                    entries.append(piece)
            return entries

        # ------------------------------------------------------------------
        # Batch over accessions
        # ------------------------------------------------------------------
        for i in range(0, len(accessions), batch_size):
            batch = accessions[i : i + batch_size]
            if not batch:
                continue

            params = {
                "format": "tsv",
                "query": " OR ".join(f"accession:{acc}" for acc in batch),
                "fields": fields,
            }

            try:
                resp = self.session.get(url, params=params, timeout=self.http_timeout)
                resp.raise_for_status()
            except Exception as e:  # pragma: no cover
                logger.error(
                    "UniProt protein annotation fetch failed for batch %d-%d: %s",
                    i + 1,
                    i + len(batch),
                    e,
                )
                continue

            lines = resp.text.splitlines()
            if len(lines) <= 1:
                continue

            for line in lines[1:]:
                if not line.strip():
                    continue
                cols = line.split("\t")
                # Pad if UniProt ever omits trailing fields
                while len(cols) < 4 + len(xref_field_list):
                    cols.append("")

                acc = cols[0].strip()
                if not acc:
                    continue

                bp_raw = cols[1]
                mf_raw = cols[2]
                cc_raw = cols[3]

                ann: Dict[str, List[str]] = {
                    "GO_BP": self._uniq(_split_go_field(bp_raw)),
                    "GO_MF": self._uniq(_split_go_field(mf_raw)),
                    "GO_CC": self._uniq(_split_go_field(cc_raw)),
                }

                # Parse all configured xref_* columns into labelled lists
                for offset, field_name in enumerate(xref_field_list):
                    raw_val = cols[4 + offset]
                    label = xref_ontology_fields[field_name]
                    ann[label] = self._uniq(_split_db_field(raw_val))

                # Only store if we have at least something
                if any(ann.values()):
                    annotations[acc] = ann

        return annotations

    def _map_refseq_to_uniprot(
        self, refseq_ids: List[str]
    ) -> Dict[str, Dict[str, List[str]]]:
        """
        Map RefSeq protein accessions to UniProtKB using bioservices.UniProt.mapping.

        Returns:
            dict:
                {refseq: {"primaryAccession": [...], "uniProtkbId": [...]}}

        NOTE:
            This uses UniProt's ID mapping job API under the hood via bioservices.UniProt.mapping.
        """
        if self.uniprot is None:
            return {}

        mapping: Dict[str, Dict[str, List[str]]] = {}

        # bioservices.UniProt.mapping can take lists directly
        try:
            res = self.uniprot.mapping(
                fr="RefSeq_Protein",
                to="UniProtKB",
                query=refseq_ids,
                max_waiting_time=600,
            )
        except Exception as e:  # pragma: no cover
            logger.error("UniProt mapping failed: %s", e)
            return {}

        results = res.get("results", []) if isinstance(res, dict) else []
        for item in results:
            src = item.get("from")
            to = item.get("to", {})
            primary = to.get("primaryAccession")
            up_id = to.get("uniProtkbId")
            if not src:
                continue
            if src not in mapping:
                mapping[src] = {"primaryAccession": [], "uniProtkbId": []}
            if primary:
                mapping[src]["primaryAccession"].append(primary)
            if up_id:
                mapping[src]["uniProtkbId"].append(up_id)

        return mapping

    def enrich_proteins(
        self,
        input_name: str = "Protein_Properties_Processed.csv",
        output_name: str = "Protein_Properties_WithOntologies.csv",
        accession_col: str = "ProteinRefSeqAccession",
    ) -> Optional[pd.DataFrame]:
        """
        Enrich Protein nodes with UniProt accession, GO terms, and a broad set
        of ontology / pathway / domain / interaction / disease cross-references.

        Input:
            - CSV with a column `accession_col` containing RefSeq or UniProt
              protein accessions.

        Output columns added:
            UniProt_PrimaryAccession
            UniProt_IDs

            GO_BP_Protein_IDs
            GO_MF_Protein_IDs
            GO_CC_Protein_IDs

            InterPro_IDs
            Pfam_IDs
            PROSITE_IDs
            PANTHER_Family_IDs
            PIRSF_Family_IDs
            PRINTS_Family_IDs
            Gene3D_Domain_IDs
            CATH_FunFam_IDs
            SFLD_Family_IDs
            SUPFAM_IDs
            SMART_Domain_IDs
            CDD_Domain_IDs
            HAMAP_Profile_IDs
            PRO_IDs

            PDB_Structure_IDs
            Reactome_Pathway_IDs
            KEGG_Pathway_IDs
            UniPathway_IDs
            PlantReactome_Pathway_IDs

            STRING_IDs
            BioGRID_IDs
            IntAct_IDs

            NCBI_Gene_IDs
            Ensembl_Protein_IDs
            Orphanet_Disease_IDs
            OMIM_IDs
        """
        input_path = f"{self.data_dir}/{input_name}"

        try:
            df = pd.read_csv(input_path, low_memory=False)
        except FileNotFoundError:
            logger.info("Protein input file not found: %s (skipping)", input_path)
            return None

        logger.info("Running protein ontology enrichment (UniProt/GO/xref ontologies)...")

        if self.uniprot is None:
            logger.warning(
                "Skipping protein UniProt/GO/ontology enrichment "
                "(bioservices.UniProt not available)."
            )
            out_path = f"{self.data_dir}/{output_name}"
            df.to_csv(out_path, index=False)
            logger.info("Saved (unenriched) protein table to %s", out_path)
            return df

        if accession_col not in df.columns:
            logger.error("Column '%s' not found in %s", accession_col, input_path)
            return None

        # ------------------------------------------------------
        # 1) Map RefSeq → UniProt (where possible)
        # ------------------------------------------------------
        raw_ids = (
            df[accession_col]
            .dropna()
            .astype(str)
            .str.strip()
            .unique()
            .tolist()
        )

        logger.info("Mapping %d unique accessions to UniProtKB", len(raw_ids))

        # If all IDs already look like UniProt accessions, skip RefSeq mapping
        if raw_ids and all(self._looks_like_uniprot_accession(a) for a in raw_ids):
            logger.info(
                "Accession column '%s' appears to contain UniProt IDs; "
                "using them directly as primary accessions.",
                accession_col,
            )
            refseq_to_uniprot = {
                acc: {"primaryAccession": [acc], "uniProtkbId": [acc]}
                for acc in raw_ids
            }
        else:
            refseq_to_uniprot = self._map_refseq_to_uniprot(raw_ids)

        # First pass: assign UniProt_PrimaryAccession / UniProt_IDs
        up_primary: List[Optional[str]] = []
        up_ids: List[str] = []
        primary_accs_set: set = set()

        for _, row in df.iterrows():
            raw_acc = (
                str(row[accession_col]).strip()
                if not pd.isna(row[accession_col])
                else ""
            )
            mapping_entry = refseq_to_uniprot.get(raw_acc, {})

            prim_raw = mapping_entry.get("primaryAccession", [])
            ids_raw = mapping_entry.get("uniProtkbId", [])

            # Normalise to lists
            if isinstance(prim_raw, str):
                prim_list = [prim_raw]
            elif isinstance(prim_raw, list):
                prim_list = prim_raw
            else:
                prim_list = []

            if isinstance(ids_raw, str):
                ids_list = [ids_raw]
            elif isinstance(ids_raw, list):
                ids_list = ids_raw
            else:
                ids_list = []

            prim_list = [p for p in prim_list if p]
            ids_list = [u for u in ids_list if u]

            primary_acc: Optional[str] = None

            if prim_list:
                primary_acc = prim_list[0]
            else:
                if raw_acc and self._looks_like_uniprot_accession(raw_acc):
                    primary_acc = raw_acc

            if not ids_list and primary_acc:
                ids_list = [primary_acc]

            up_primary.append(primary_acc)
            if primary_acc:
                primary_accs_set.add(primary_acc)

            up_ids.append("|".join(sorted(set(ids_list))) if ids_list else "")

        df["UniProt_PrimaryAccession"] = up_primary
        df["UniProt_IDs"] = up_ids

        # ------------------------------------------------------
        # 2) Fetch GO / ontology xrefs for all primary accessions
        # ------------------------------------------------------
        primary_accs = sorted(primary_accs_set)
        logger.info(
            "Fetching UniProt annotations for %d proteins",
            len(primary_accs),
        )

        prot_ann = self._fetch_uniprot_protein_annotations(primary_accs)

        # ------------------------------------------------------
        # 3) Second pass: attach ontology columns using primary accession
        # ------------------------------------------------------
        go_bp_ids: List[str] = []
        go_mf_ids: List[str] = []
        go_cc_ids: List[str] = []

        # Map annotation labels -> output column names
        xref_output_map: Dict[str, str] = {
            "InterPro": "InterPro_IDs",
            "Pfam": "Pfam_IDs",
            "PROSITE": "PROSITE_IDs",
            "PANTHER": "PANTHER_Family_IDs",
            "PIRSF": "PIRSF_Family_IDs",
            "PRINTS": "PRINTS_Family_IDs",
            "Gene3D": "Gene3D_Domain_IDs",
            "FunFam": "CATH_FunFam_IDs",
            "SFLD": "SFLD_Family_IDs",
            "SUPFAM": "SUPFAM_IDs",
            "SMART": "SMART_Domain_IDs",
            "CDD": "CDD_Domain_IDs",
            "HAMAP": "HAMAP_Profile_IDs",
            "PRO": "PRO_IDs",

            "PDB": "PDB_Structure_IDs",
            "Reactome": "Reactome_Pathway_IDs",
            "KEGG": "KEGG_Pathway_IDs",
            "UniPathway": "UniPathway_IDs",
            "PlantReactome": "PlantReactome_Pathway_IDs",

            "STRING": "STRING_IDs",
            "BioGRID": "BioGRID_IDs",
            "IntAct": "IntAct_IDs",

            "NCBIGene": "NCBI_Gene_IDs",
            "Ensembl": "Ensembl_Protein_IDs",
            "Orphanet": "Orphanet_Disease_IDs",
            "MIM": "OMIM_IDs",
        }

        # Prepare lists to collect values per column
        xref_values_lists: Dict[str, List[str]] = {
            col_name: [] for col_name in xref_output_map.values()
        }

        def _join(xs: List[str]) -> str:
            xs = [x for x in xs if x]
            return "|".join(sorted(set(xs))) if xs else ""

        for _, row in df.iterrows():
            acc = (
                str(row["UniProt_PrimaryAccession"]).strip()
                if not pd.isna(row["UniProt_PrimaryAccession"])
                else ""
            )
            ann = prot_ann.get(acc, {}) if acc else {}

            go_bp_ids.append(_join(ann.get("GO_BP", [])))
            go_mf_ids.append(_join(ann.get("GO_MF", [])))
            go_cc_ids.append(_join(ann.get("GO_CC", [])))

            for label, col_name in xref_output_map.items():
                xref_values_lists[col_name].append(_join(ann.get(label, [])))

        df["GO_BP_Protein_IDs"] = go_bp_ids
        df["GO_MF_Protein_IDs"] = go_mf_ids
        df["GO_CC_Protein_IDs"] = go_cc_ids

        for col_name, values in xref_values_lists.items():
            df[col_name] = values

        # ------------------------------------------------------
        # 4) Ontology counts summary + transparency metadata
        # ------------------------------------------------------
        ontology_cols = [
            "UniProt_PrimaryAccession",
            "UniProt_IDs",
            "GO_BP_Protein_IDs",
            "GO_MF_Protein_IDs",
            "GO_CC_Protein_IDs",
        ] + list(xref_values_lists.keys())

        ontology_counts: Dict[str, int] = {}
        for col in ontology_cols:
            ontology_counts[col] = self._count_non_empty_column(df[col])

        num_rows = len(df)
        ontology_percent = {
            k: (v / num_rows * 100.0 if num_rows > 0 else 0.0)
            for k, v in ontology_counts.items()
        }

        # Rows with no ontology at all (except possibly UniProt_IDs)
        ontology_presence_cols = [
            "UniProt_PrimaryAccession",
            "GO_BP_Protein_IDs",
            "GO_MF_Protein_IDs",
            "GO_CC_Protein_IDs",
        ] + [
            c for c in xref_values_lists.keys()
        ]

        missing_all: List[str] = []
        for _, row in df.iterrows():
            has_any = False
            for col in ontology_presence_cols:
                val = row.get(col)
                if pd.notna(val) and str(val).strip() != "":
                    has_any = True
                    break

            if not has_any:
                acc_val = (
                    str(row[accession_col]).strip()
                    if not pd.isna(row[accession_col])
                    else ""
                )
                if acc_val:
                    missing_all.append(acc_val)

        num_missing = len(missing_all)
        num_success = num_rows - num_missing

        logger.info("=== Protein Ontology Retrieval Summary ===")
        logger.info("Total proteins processed: %d", num_rows)
        logger.info("Successful hits (any ontology): %d", num_success)
        logger.info("Missing entries (no ontologies at all): %d", num_missing)
        logger.info("Ontology counts / coverage:")
        for k, v in ontology_counts.items():
            logger.info("  %s: %d  (%.2f%%)", k, v, ontology_percent[k])

        run_summary = {
            "total_proteins_processed": int(num_rows),
            "total_unique_input_accessions": int(len(raw_ids)),
            "total_successful_hits": int(num_success),
            "total_missing_entries": int(num_missing),
            "missing_input_accessions": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "retrieved_fields": ontology_cols,
            "input_file": input_path,
            "output_file": f"{self.data_dir}/{output_name}",
            "accession_column": accession_col,
            "uniprot_fields": [
                "accession",
                "go_p",
                "go_f",
                "go_c",
            ] + list(prot_ann[next(iter(prot_ann))].keys() - {"GO_BP", "GO_MF", "GO_CC"})
              if prot_ann
              else [],
            "created_at": datetime.datetime.utcnow().isoformat() + "Z",
            "batch_size_uniprot": 200,
            "source": "UniProtKB REST",
        }

        summary_path = f"{self.data_dir}/Protein_Ontology_RunSummary.json"
        self._save_run_summary(run_summary, summary_path)

        # ------------------------------------------------------
        # 5) Save enriched version
        # ------------------------------------------------------
        df["OntologyEnrichedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.ProteinUniProt/2.0"
        df["OntologySources"] = "UniProtKB"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved enriched protein table to %s", out_path)

        return df

    # ------------------------------------------------------------------
    # Keyword extraction from assay descriptions
    # ------------------------------------------------------------------
    def _init_nlp(self) -> None:
        """
        Lazy-load spaCy English model if available.
        This gives us POS tags, lemmas, and built-in stopword handling
        without hard-coding any word lists here.
        """
        if self._nlp is not None:
            return

        try:
            import spacy
            # you will need: python -m spacy download en_core_web_sm
            self._nlp = spacy.load("en_core_web_sm")
            logger.info("Loaded spaCy model 'en_core_web_sm' for keyword extraction.")
        except Exception as e:  # pragma: no cover
            logger.warning(
                "spaCy or 'en_core_web_sm' not available (%s). "
                "Falling back to simple heuristic keyword extraction.", e
            )
            self._nlp = None

    @staticmethod
    def _extract_keywords_fallback(text: str, max_keywords: int = 10) -> Set[str]:
        """
        Very simple, fully unsupervised keyword extractor.

        - No stopword list at all.
        - Tokenise on alphanumerics.
        - Drop very short or purely numeric tokens.
        - Score by frequency and length.

        This is *not* truly semantic, but provides a sensible ranking
        without manually listing words.
        """
        if not text:
            return set()

        tokens = re.findall(r"[A-Za-z0-9]+", text.lower())
        if not tokens:
            return set()

        filtered: List[str] = []
        for t in tokens:
            if t.isdigit():
                continue
            if len(t) <= 2:
                continue
            filtered.append(t)

        if not filtered:
            return set()

        freq = Counter(filtered)

        scores: Dict[str, float] = {}
        for tok, count in freq.items():
            length = len(tok)
            # favour slightly longer tokens, but cap to avoid crazy bias
            length_factor = min(max((length - 3) / 7.0, 0.0), 1.0)
            scores[tok] = float(count) * (1.0 + length_factor)

        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        top_tokens = [tok for tok, _ in ranked[:max_keywords]]

        return set(top_tokens)

    def _extract_keywords(self, text: str, max_keywords: int = 10) -> Set[str]:
        """
        Keyword extractor that aims for semantic importance:

        Preferred path (if spaCy is available):
            - Run spaCy pipeline.
            - Keep content-bearing lemmas (NOUN, PROPN, ADJ, VERB).
            - Rely on spaCy's internal linguistic knowledge (including
              its own stopword handling) instead of any manual word list.
            - Score by frequency + positional weight.

        Fallback path (no spaCy):
            - Call `_extract_keywords_fallback` (no stopword list).
        """
        if not text:
            return set()

        # Try to use spaCy if available
        if getattr(self, "_nlp", None) is None:
            self._init_nlp()

        if self._nlp is not None:
            doc = self._nlp(text)

            candidates: List[tuple[str, float]] = []
            doc_len = max(len(doc), 1)

            for i, token in enumerate(doc):
                # Use spaCy's own idea of "stop words" and punctuation,
                # but we do NOT hard-code any list in our own code.
                if token.is_stop or token.is_punct or token.is_space:
                    continue

                # Focus on content-bearing parts of speech
                if token.pos_ not in ("NOUN", "PROPN", "ADJ", "VERB"):
                    continue

                lemma = token.lemma_.lower().strip()
                if not lemma:
                    continue
                if lemma.isdigit():
                    continue
                if len(lemma) <= 2:
                    continue

                # Simple positional weighting: earlier tokens weigh slightly more
                pos_weight = 1.0 + 0.5 * (1.0 - i / doc_len)
                candidates.append((lemma, pos_weight))

            if candidates:
                scores: Counter = Counter()
                for lemma, w in candidates:
                    scores[lemma] += w

                top_lemmas = [tok for tok, _ in scores.most_common(max_keywords)]
                return set(top_lemmas)

            # If spaCy finds nothing, fall back to unsupervised variant
            return self._extract_keywords_fallback(text, max_keywords)

        # No spaCy: unsupervised fallback
        return self._extract_keywords_fallback(text, max_keywords)

    @staticmethod
    def _build_search_text(row: pd.Series) -> str:
        """
        Build a compact text for ontology lookup from the most informative
        columns (short and descriptive):

        Priority:
            AssayActivityName > AssayType > AssayName > AssayDescription

        Also:
        - Strip generic 'Other'
        - Truncate to 120 chars to avoid confusing OLS with very long blobs
        """
        candidates = []
        for col in ["AssayActivityName", "AssayType", "AssayName"]:
            if col in row and pd.notna(row[col]):
                txt = str(row[col]).strip()
                if txt:
                    candidates.append(txt)

        if not candidates and "AssayDescription" in row and pd.notna(row["AssayDescription"]):
            candidates.append(str(row["AssayDescription"]).strip())

        if not candidates:
            return ""

        text = " ; ".join(candidates)

        # Very lightweight cleaning
        text = text.replace("Other ;", "").replace("other ;", "").strip()

        # Keep it reasonably short for OLS
        return text[:120]

    def enrich_assays(
        self,
        input_name: str = "Assay_Properties_Processed.csv",
        output_name: str = "Assay_Properties_WithOntologies.csv",
        name_cols: Optional[List[str]] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Enrich BioAssay nodes with assay ontologies via OLS.

        Ontologies:
            - BAO (BioAssay Ontology)  -> primary assay type
            - OBI (Ontology for Biomedical Investigations) -> method-level context

        Strategy:
            - Build a 'search string' per assay from a combination of columns
            - Query OLS in each ontology (bao, obi)
            - For each ontology:
                * pick best match by string_score
                * pick best match by keyword_score
                * pick best match by combined_score
            - Store all three “views” in separate columns
            - Keep legacy BAO_ID / BAO_Label / BAO_Score as aliases of the
            combined-score mapping for backwards compatibility.

        Output columns added (per ontology = BAO, OBI):

            *_ID_String, *_IRI_String, *_Label_String, *_StringScore
            *_ID_Keyword, *_IRI_Keyword, *_Label_Keyword, *_KeywordScore
            *_ID, *_IRI, *_Label, *_Score   (combined / main mapping)

        A JSON summary is written to:
            Data/Nodes/Assay_Ontology_RunSummary.json
        """
        input_path = f"{self.data_dir}/{input_name}"

        try:
            df = pd.read_csv(input_path, low_memory=False)
            # df = df.iloc[:, :3]
        except FileNotFoundError:
            logger.info("Assay input file not found: %s (skipping)", input_path)
            return None

        if name_cols is None:
            name_cols = ["AssayName", "AssayActivityName", "AssayType", "AssayDescription"]

        # ------------------------------------------------------------------
        # Build search text per row from *informative* name/description cols
        # ------------------------------------------------------------------
        search_texts = [self._build_search_text(row) for _, row in df.iterrows()]

        # Helper to convert an OBO-style IRI to CURIE-like ID (e.g. BAO_0000015 -> BAO:0000015)
        def _iri_to_curie(iri: str) -> str:
            if not iri:
                return ""
            local = iri.rsplit("/", 1)[-1]
            # If a fragment is used
            if "#" in local:
                local = local.split("#")[-1]
            if "_" in local:
                prefix, rest = local.split("_", 1)
                return f"{prefix}:{rest}"
            return local

        # ------------------------------------------------------------------
        # Cache OLS lookups per (ontology, text) to avoid redundant calls
        # Each cache entry: {
        #   "string":   (iri, label, score),
        #   "keyword":  (iri, label, score),
        #   "combined": (iri, label, score),
        # }
        # ------------------------------------------------------------------
        target_ontologies = ["bao", "obi"]
        cache: Dict[str, Dict[str, Dict[str, Tuple[str, str, float]]]] = {
            ont: {} for ont in target_ontologies
        }

        def _lookup_ontology_multi(text: str, ontology: str) -> Dict[str, Tuple[str, str, float]]:
            """
            For a given (text, ontology), compute and return the *best* match
            according to three different scoring strategies:

                - "string":   max string_score
                - "keyword":  max keyword_score
                - "combined": max combined_score (subject to threshold)

            Returns a dict:
                {
                "string":   (iri, label, string_score),
                "keyword":  (iri, label, keyword_score),
                "combined": (iri, label, combined_score),
                }

            Empty/failed cases return empty IRIs/labels and 0.0 scores.
            """
            if not text:
                return {
                    "string":   ("", "", 0.0),
                    "keyword":  ("", "", 0.0),
                    "combined": ("", "", 0.0),
                }

            cache_for_ont = cache[ontology]
            if text in cache_for_ont:
                return cache_for_ont[text]

            term_clean = str(text).strip()
            if not term_clean:
                result = {
                    "string":   ("", "", 0.0),
                    "keyword":  ("", "", 0.0),
                    "combined": ("", "", 0.0),
                }
                cache_for_ont[text] = result
                return result

            url = f"{self.ols_base_url}/search"
            params = {
                "q": term_clean,
                "ontology": ontology,
                "rows": 10,  # a few more candidates to let scores differentiate
                "queryFields": "label,synonym,short_form",
            }

            data = None
            for attempt in range(1, 4):
                try:
                    resp = self.session.get(url, params=params, timeout=self.http_timeout)
                    resp.raise_for_status()
                    data = resp.json()
                    break
                except Exception as e:  # pragma: no cover
                    logger.warning(
                        "OLS query failed for '%s' in %s (attempt %d/3): %s",
                        term_clean,
                        ontology,
                        attempt,
                        e,
                    )
                    if attempt == 3:
                        logger.error(
                            "OLS query permanently failed for '%s' in %s after %d attempts",
                            term_clean,
                            ontology,
                            attempt,
                        )
                        result = {
                            "string":   ("", "", 0.0),
                            "keyword":  ("", "", 0.0),
                            "combined": ("", "", 0.0),
                        }
                        cache_for_ont[text] = result
                        return result
                    time.sleep(0.5 * attempt)

            docs = (data or {}).get("response", {}).get("docs", [])
            if not docs:
                logger.debug(
                    "OLS returned no docs for '%s' in %s (raw data keys: %s)",
                    term_clean, ontology, list((data or {}).keys())
                )
                result = {
                    "string":   ("", "", 0.0),
                    "keyword":  ("", "", 0.0),
                    "combined": ("", "", 0.0),
                }
                cache_for_ont[text] = result
                return result

            term_lower = term_clean.lower()

            # Pre-compute query keywords once
            src_keywords = self._extract_keywords(term_clean)

            # Best candidates per scoring strategy
            best_string = ("", "", 0.0)    # (iri, label, string_score)
            best_keyword = ("", "", 0.0)   # (iri, label, keyword_score)
            best_combined = ("", "", 0.0)  # (iri, label, combined_score)

            for d in docs:
                label = (d.get("label") or term_clean).strip()
                label_lower = label.lower()

                synonyms = d.get("synonym") or []
                if isinstance(synonyms, str):
                    synonyms = [synonyms]

                short_form = d.get("short_form") or []
                if isinstance(short_form, str):
                    short_form = [short_form]

                # -----------------------------
                # 1) String similarity
                # -----------------------------
                # exact match boost
                if label_lower == term_lower or any(
                    str(s).strip().lower() == term_lower for s in synonyms
                ):
                    string_score = 1.0
                else:
                    string_score = SequenceMatcher(
                        None,
                        term_lower,
                        label_lower
                    ).ratio()

                # -----------------------------
                # 2) Keyword overlap
                # -----------------------------
                if src_keywords:
                    target_text = " ".join([label] + synonyms + short_form)
                    tgt_keywords = self._extract_keywords(target_text)
                    if tgt_keywords:
                        intersect = src_keywords & tgt_keywords
                        keyword_score = len(intersect) / len(src_keywords)
                    else:
                        keyword_score = 0.0
                else:
                    keyword_score = 0.0

                # -----------------------------
                # 3) Combined score
                # -----------------------------
                alpha = 0.7  # favour string a bit
                combined_score = alpha * string_score + (1.0 - alpha) * keyword_score

                iri = d.get("iri") or ""

                # Update best-by-string
                if string_score > best_string[2]:
                    best_string = (iri, label, float(string_score))

                # Update best-by-keyword
                if keyword_score > best_keyword[2]:
                    best_keyword = (iri, label, float(keyword_score))

                # Update best-by-combined
                if combined_score > best_combined[2]:
                    best_combined = (iri, label, float(combined_score))

            # Apply ontology-specific minimum threshold on the combined mapping only
            default_min_scores = {"bao": 0.25, "obi": 0.30}
            ontology_min_scores = getattr(self, "ontology_min_scores", default_min_scores)
            eff_min_score = float(ontology_min_scores.get(ontology, 0.1))

            iri_c, label_c, combined_score = best_combined
            if combined_score < eff_min_score:
                # If combined score is too low, blank out the main mapping.
                best_combined = ("", "", 0.0)

            result = {
                "string":   best_string,
                "keyword":  best_keyword,
                "combined": best_combined,
            }
            cache_for_ont[text] = result

            # Respect rate limiting if configured
            if self.OLS_SLEEP_SECONDS and self.OLS_SLEEP_SECONDS > 0:
                time.sleep(self.OLS_SLEEP_SECONDS)

            return result

        # ------------------------------------------------------------------
        # Perform enrichment
        # ------------------------------------------------------------------
        # BAO columns
        bao_id_string: List[str] = []
        bao_iri_string: List[str] = []
        bao_label_string: List[str] = []
        bao_string_score: List[float] = []

        bao_id_keyword: List[str] = []
        bao_iri_keyword: List[str] = []
        bao_label_keyword: List[str] = []
        bao_keyword_score: List[float] = []

        bao_ids: List[str] = []         # combined / main
        bao_iris: List[str] = []
        bao_labels: List[str] = []
        bao_scores: List[float] = []    # combined_score

        # OBI columns
        obi_id_string: List[str] = []
        obi_iri_string: List[str] = []
        obi_label_string: List[str] = []
        obi_string_score: List[float] = []

        obi_id_keyword: List[str] = []
        obi_iri_keyword: List[str] = []
        obi_label_keyword: List[str] = []
        obi_keyword_score: List[float] = []

        obi_ids: List[str] = []         # combined / main
        obi_iris: List[str] = []
        obi_labels: List[str] = []
        obi_scores: List[float] = []    # combined_score

        logger.info("Enriching %d assays with BAO/OBI via OLS (per-score mappings)", len(df))

        for text in search_texts:
            # ---------------- BAO ----------------
            bao_res = _lookup_ontology_multi(text, "bao")

            # best by string
            b_iri_s, b_label_s, b_score_s = bao_res["string"]
            bao_iri_string.append(b_iri_s)
            bao_label_string.append(b_label_s)
            bao_id_string.append(_iri_to_curie(b_iri_s))
            bao_string_score.append(b_score_s)

            # best by keyword
            b_iri_k, b_label_k, b_score_k = bao_res["keyword"]
            bao_iri_keyword.append(b_iri_k)
            bao_label_keyword.append(b_label_k)
            bao_id_keyword.append(_iri_to_curie(b_iri_k))
            bao_keyword_score.append(b_score_k)

            # best by combined (main mapping)
            b_iri_c, b_label_c, b_score_c = bao_res["combined"]
            bao_iris.append(b_iri_c)
            bao_labels.append(b_label_c)
            bao_ids.append(_iri_to_curie(b_iri_c))
            bao_scores.append(b_score_c)

            # ---------------- OBI ----------------
            obi_res = _lookup_ontology_multi(text, "obi")

            # best by string
            o_iri_s, o_label_s, o_score_s = obi_res["string"]
            obi_iri_string.append(o_iri_s)
            obi_label_string.append(o_label_s)
            obi_id_string.append(_iri_to_curie(o_iri_s))
            obi_string_score.append(o_score_s)

            # best by keyword
            o_iri_k, o_label_k, o_score_k = obi_res["keyword"]
            obi_iri_keyword.append(o_iri_k)
            obi_label_keyword.append(o_label_k)
            obi_id_keyword.append(_iri_to_curie(o_iri_k))
            obi_keyword_score.append(o_score_k)

            # best by combined (main mapping)
            o_iri_c, o_label_c, o_score_c = obi_res["combined"]
            obi_iris.append(o_iri_c)
            obi_labels.append(o_label_c)
            obi_ids.append(_iri_to_curie(o_iri_c))
            obi_scores.append(o_score_c)

        # ------------------------------------------------------------------
        # Attach ontology columns
        #   - String-based mapping
        #   - Keyword-based mapping
        #   - Combined (main) mapping (keeps legacy column names)
        # ------------------------------------------------------------------
        # BAO
        df["BAO_ID_String"] = bao_id_string
        df["BAO_IRI_String"] = bao_iri_string
        df["BAO_Label_String"] = bao_label_string
        df["BAO_StringScore"] = bao_string_score

        df["BAO_ID_Keyword"] = bao_id_keyword
        df["BAO_IRI_Keyword"] = bao_iri_keyword
        df["BAO_Label_Keyword"] = bao_label_keyword
        df["BAO_KeywordScore"] = bao_keyword_score

        # Legacy / main: combined mapping
        df["BAO_ID"] = bao_ids
        df["BAO_IRI"] = bao_iris
        df["BAO_Label"] = bao_labels
        df["BAO_Score"] = bao_scores

        # OBI
        df["OBI_ID_String"] = obi_id_string
        df["OBI_IRI_String"] = obi_iri_string
        df["OBI_Label_String"] = obi_label_string
        df["OBI_StringScore"] = obi_string_score

        df["OBI_ID_Keyword"] = obi_id_keyword
        df["OBI_IRI_Keyword"] = obi_iri_keyword
        df["OBI_Label_Keyword"] = obi_label_keyword
        df["OBI_KeywordScore"] = obi_keyword_score

        # Legacy / main: combined mapping
        df["OBI_ID"] = obi_ids
        df["OBI_IRI"] = obi_iris
        df["OBI_Label"] = obi_labels
        df["OBI_Score"] = obi_scores

        # ------------------------------------------------------------------
        # Ontology counts summary + transparency metadata
        # (based on main / combined mapping)
        # ------------------------------------------------------------------
        ontology_cols = ["BAO_ID", "OBI_ID"]

        ontology_counts: Dict[str, int] = {}
        for col in ontology_cols:
            ontology_counts[col] = self._count_non_empty_column(df[col])

        num_rows = len(df)
        ontology_percent = {
            k: (v / num_rows * 100.0 if num_rows > 0 else 0.0)
            for k, v in ontology_counts.items()
        }

        # Identify assays with no ontology at all (BAO and OBI both empty for main mapping)
        missing_all: List[int] = []
        id_col = "AssayID" if "AssayID" in df.columns else None

        for _, row in df.iterrows():
            has_any = any(
                pd.notna(row[col]) and str(row[col]).strip() != ""
                for col in ontology_cols
            )
            if not has_any and id_col is not None:
                aid = row.get(id_col)
                if pd.notna(aid):
                    missing_all.append(int(aid))

        num_missing = len(missing_all)
        num_success = num_rows - num_missing

        logger.info("=== Assay Ontology Retrieval Summary ===")
        logger.info("Total assays processed: %d", num_rows)
        logger.info("Successful hits (any of BAO/OBI, combined score): %d", num_success)
        logger.info("Missing entries (no BAO/OBI at all): %d", num_missing)
        logger.info("Ontology counts / coverage (by ID, combined mapping):")
        for k, v in ontology_counts.items():
            logger.info("  %s: %d  (%.2f%%)", k, v, ontology_percent[k])

        # Build and save JSON summary
        run_summary = {
            "total_assays_processed": int(num_rows),
            "total_successful_hits_any": int(num_success),
            "total_missing_entries": int(num_missing),
            "missing_assay_ids": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "retrieved_fields": [
                # BAO
                "BAO_ID_String", "BAO_IRI_String", "BAO_Label_String", "BAO_StringScore",
                "BAO_ID_Keyword", "BAO_IRI_Keyword", "BAO_Label_Keyword", "BAO_KeywordScore",
                "BAO_ID", "BAO_IRI", "BAO_Label", "BAO_Score",
                # OBI
                "OBI_ID_String", "OBI_IRI_String", "OBI_Label_String", "OBI_StringScore",
                "OBI_ID_Keyword", "OBI_IRI_Keyword", "OBI_Label_Keyword", "OBI_KeywordScore",
                "OBI_ID", "OBI_IRI", "OBI_Label", "OBI_Score",
            ],
            "input_file": input_path,
            "output_file": f"{self.data_dir}/{output_name}",
            "name_columns_used": name_cols,
            "ols_ontologies": ["bao", "obi"],
            "ols_base_url": self.ols_base_url,
            "created_at": datetime.datetime.utcnow().isoformat() + "Z",
            "ols_sleep_seconds": self.OLS_SLEEP_SECONDS,
            "source": "EBI OLS (BAO;OBI)",
        }

        summary_path = f"{self.data_dir}/Assay_Ontology_RunSummary.json"
        self._save_run_summary(run_summary, summary_path)

        # ------------------------------------------------------------------
        # Save enriched version
        # ------------------------------------------------------------------
        enrichment_time = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnrichedAt"] = enrichment_time
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.AssaysOLS/1.3"
        df["OntologySources"] = "OLS(BAO;OBI)"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved assay-enriched table to %s", out_path)

        return df
