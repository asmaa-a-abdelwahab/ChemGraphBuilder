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

from collections import Counter, deque
import datetime
import logging
import re
import time
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union
from biothings_client import get_client
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from difflib import SequenceMatcher
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

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

    # External ontology endpoints used for compound enrichment
    PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov"
    MESH_BASE_URL = "https://id.nlm.nih.gov/mesh"
    CHEBI_BASE_URL = "https://www.ebi.ac.uk/chebi/backend/api/public"

    # Token extractors for IDs embedded in PubChem PUG-View records
    _MESH_ID_TOKEN_RE = re.compile(r"\b([DCM]\d{6,9})\b")
    _CHEBI_TOKEN_RE = re.compile(r"\bCHEBI:\d+\b")

    def __init__(
        self,
        data_dir: str = "Data/Nodes",
        species: str = "human",
        mygene_email: Optional[str] = None,
        ols_base_url: str = "https://www.ebi.ac.uk/ols4/api",
        http_timeout: int = 30,
        # PubChem tuning (avoid 504s + control run time)
        pubchem_xrefs_batch_size: int = 50,
        pubchem_xrefs_min_batch_size: int = 10,
        pubchem_xrefs_max_retries: int = 3,
        pubchem_xrefs_retry_sleep: float = 5.0,
        pubchem_polite_sleep: float = 0.2,
    ):
        self._nlp = None   # lazy-initialised spaCy model
        self.data_dir = data_dir.rstrip("/")
        self.species = species
        self.ols_base_url = ols_base_url
        self.http_timeout = http_timeout

        # PubChem runtime controls (used by batch xrefs + polite sleeps)
        self.pubchem_xrefs_batch_size = int(pubchem_xrefs_batch_size)
        self.pubchem_xrefs_min_batch_size = int(pubchem_xrefs_min_batch_size)
        self.pubchem_xrefs_max_retries = int(pubchem_xrefs_max_retries)
        self.pubchem_xrefs_retry_sleep = float(pubchem_xrefs_retry_sleep)
        self.pubchem_polite_sleep = float(pubchem_polite_sleep)

        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "ChemGraphBuilder-NodesOntologyEnricher/1.0"})

        # Connection pooling + automatic retries (faster + more resilient for large runs)
        retry_cfg = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=0.5,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=("GET", "POST"),
            raise_on_status=False,
        )
        adapter = HTTPAdapter(max_retries=retry_cfg, pool_connections=64, pool_maxsize=64)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)

        # Internal caches (in-memory) to minimize repeated HTTP calls
        self._pug_view_cache: Dict[int, Optional[Dict[str, Any]]] = {}
        self._chebi_parents_cache: Dict[str, List[Dict[str, str]]] = {}

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
        # Optional FigureB auxiliary nodes (auto-skip if input files don't exist)
        self.enrich_experimental_contexts()
        self.enrich_assay_endpoints()


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


    # ------------------------------------------------------------------
    # Compound enrichment helpers: PubChem PUG-View, MeSH lookup, ChEBI parents
    # ------------------------------------------------------------------
    def _safe_get_json(self, url: str, params: Optional[dict] = None, timeout: Optional[int] = None) -> Optional[Any]:
        """Best-effort JSON GET that returns None on HTTP errors/timeouts."""
        try:
            r = self.session.get(
                url,
                params=params,
                timeout=timeout or self.http_timeout,
                headers={"Accept": "application/json"},
            )
            # PubChem may return 404 for missing sections; treat as empty.
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except Exception:
            return None

    def _pubchem_pug_view_record(self, cid: int) -> Optional[Dict[str, Any]]:
        """Fetch full PubChem PUG-View record for a CID (cached)."""
        try:
            cid = int(cid)
        except Exception:
            return None

        if cid in self._pug_view_cache:
            return self._pug_view_cache[cid]

        url = f"{self.PUBCHEM_BASE_URL}/rest/pug_view/data/compound/{cid}/JSON"
        data = self._safe_get_json(url, params={"response_type": "display"})
        rec = None
        if isinstance(data, dict):
            rec_obj = data.get("Record")
            if isinstance(rec_obj, dict):
                rec = rec_obj

        self._pug_view_cache[cid] = rec
        return rec

    @staticmethod
    def _walk_sections(section_list: Iterable[Dict[str, Any]]) -> Iterable[Dict[str, Any]]:
        """Depth-first walk over PubChem PUG-View Record.Section nodes."""
        stack = list(section_list)[::-1]
        while stack:
            sec = stack.pop()
            yield sec
            for child in (sec.get("Section") or [])[::-1]:
                if isinstance(child, dict):
                    stack.append(child)

    @staticmethod
    def _extract_strings_from_value(value_obj: Any) -> List[str]:
        """Extract plain strings from PubChem PUG-View Value objects."""
        out: List[str] = []
        if isinstance(value_obj, dict):
            swm = value_obj.get("StringWithMarkup")
            if isinstance(swm, list):
                for item in swm:
                    if isinstance(item, dict):
                        s = item.get("String")
                        if isinstance(s, str) and s.strip():
                            out.append(s.strip())
            # Sometimes nested structures exist; be permissive
            for v in value_obj.values():
                if isinstance(v, (dict, list)) and v is not swm:
                    out.extend(NodesOntologyEnricher._extract_strings_from_value(v))
        elif isinstance(value_obj, list):
            for v in value_obj:
                out.extend(NodesOntologyEnricher._extract_strings_from_value(v))
        return out

    @classmethod
    def _extract_terms_by_toc_heading(cls, record: Dict[str, Any], toc_heading: str) -> List[str]:
        terms: Set[str] = set()
        for sec in cls._walk_sections(record.get("Section", []) or []):
            if sec.get("TOCHeading") == toc_heading:
                for info in sec.get("Information", []) or []:
                    if isinstance(info, dict):
                        terms.update(cls._extract_strings_from_value(info.get("Value")))
        return sorted(terms)

    @classmethod
    def _find_sections_heading_contains(cls, record: Dict[str, Any], substrings: Tuple[str, ...]) -> List[Dict[str, Any]]:
        hits: List[Dict[str, Any]] = []
        for sec in cls._walk_sections(record.get("Section", []) or []):
            h = sec.get("TOCHeading")
            if isinstance(h, str):
                hl = h.lower()
                if any(sub.lower() in hl for sub in substrings):
                    hits.append(sec)
        return hits

    @classmethod
    def _extract_kv_from_sections(cls, sections: List[Dict[str, Any]]) -> Dict[str, List[str]]:
        out: Dict[str, List[str]] = {}
        for sec in sections:
            for info in sec.get("Information", []) or []:
                if not isinstance(info, dict):
                    continue
                key = info.get("Name")
                if not isinstance(key, str) or not key.strip():
                    continue
                vals = cls._extract_strings_from_value(info.get("Value"))
                if not vals:
                    continue
                out.setdefault(key.strip(), [])
                for v in vals:
                    if v not in out[key.strip()]:
                        out[key.strip()].append(v)
        return out

    @staticmethod
    def _iter_all_strings(obj: Any) -> Iterable[str]:
        if isinstance(obj, dict):
            for v in obj.values():
                yield from NodesOntologyEnricher._iter_all_strings(v)
        elif isinstance(obj, list):
            for v in obj:
                yield from NodesOntologyEnricher._iter_all_strings(v)
        elif isinstance(obj, str):
            yield obj

    def _extract_ids_from_record_text(self, record: Dict[str, Any]) -> Dict[str, List[str]]:
        """Extract MeSH and ChEBI tokens embedded anywhere in a PUG-View record."""
        mesh_ids: Set[str] = set()
        chebi_ids: Set[str] = set()

        for s in self._iter_all_strings(record):
            for m in self._MESH_ID_TOKEN_RE.findall(s):
                mesh_ids.add(m)
            for c in self._CHEBI_TOKEN_RE.findall(s):
                chebi_ids.add(c)

        return {
            "mesh_unique_ids": sorted(mesh_ids),   # D*, C*, M* (best-effort)
            "chebi_ids": sorted(chebi_ids),
        }

    @staticmethod
    def _extract_pubchem_mesh_references(record: Dict[str, Any]) -> List[Dict[str, str]]:
        """Extract MeSH-related references from PubChem PUG-View record.Reference."""
        refs = record.get("Reference", []) or []
        out: List[Dict[str, str]] = []
        for ref in refs:
            if not isinstance(ref, dict):
                continue
            src = ref.get("SourceName")
            url = ref.get("URL") or ""
            if src == "Medical Subject Headings (MeSH)" or (isinstance(url, str) and "id.nlm.nih.gov/mesh/" in url):
                out.append(
                    {
                        "mesh_source_id": str(ref.get("SourceID") or ""),
                        "label": str(ref.get("Name") or ""),
                        "url": str(url or ""),
                    }
                )

        # de-dup
        seen = set()
        uniq: List[Dict[str, str]] = []
        for x in out:
            k = (x["mesh_source_id"], x["label"], x["url"])
            if k not in seen:
                seen.add(k)
                uniq.append(x)
        return uniq

    def _chebi_get_parent_nodes(self, chebi_id: str) -> List[Dict[str, str]]:
        """Fetch ChEBI ontology parents for a CHEBI ID (cached, best-effort)."""
        chebi_id = (chebi_id or "").strip()
        if not chebi_id:
            return []
        if chebi_id in self._chebi_parents_cache:
            return self._chebi_parents_cache[chebi_id]

        chebi_num = chebi_id.split(":")[-1]

        def _fetch(token: str) -> Optional[Any]:
            return self._safe_get_json(f"{self.CHEBI_BASE_URL}/ontology/parents/{token}/")

        payload = _fetch(chebi_id) or _fetch(chebi_num)
        out: List[Dict[str, str]] = []

        def rec(x: Any) -> None:
            if isinstance(x, dict):
                if ("chebiId" in x or "id" in x) and ("chebiName" in x or "name" in x):
                    cid = x.get("chebiId", x.get("id"))
                    name = x.get("chebiName", x.get("name"))
                    rel = x.get("relationship", x.get("relation", x.get("type", "")))
                    if isinstance(cid, str) and isinstance(name, str):
                        out.append({"chebiId": cid, "chebiName": name, "relationship": str(rel or "")})
                for v in x.values():
                    rec(v)
            elif isinstance(x, list):
                for v in x:
                    rec(v)

        if payload is not None:
            rec(payload)

        seen = set()
        uniq: List[Dict[str, str]] = []
        for r in out:
            k = (r.get("chebiId"), r.get("chebiName"), r.get("relationship"))
            if k not in seen:
                seen.add(k)
                uniq.append(r)

        self._chebi_parents_cache[chebi_id] = uniq
        return uniq

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
                out["MeSH_IDs"]       = self._uniq(gx("mesh"))
                out["KEGG_IDs"]     = self._uniq(gx("kegg"))
                out["HMDB_IDs"]     = self._uniq(gx("hmdb"))
                out["DTXSID_IDs"]   = self._uniq(gx("dtxsid"))
                out["DBURL"]            = self._uniq(gx("url"))
                out["CAS_RN"]    += self._uniq(gx("cas"))  # add CAS from pubchem.xrefs

                results[ik] = out

        return results

    def _fetch_compound_xrefs_pubchem_batch(
        self,
        cids: List[int],
        batch_size: int = 50,
        max_retries: int = 3,
        retry_sleep: float = 5.0,
        min_batch_size: int = 10,
    ) -> Dict[int, Dict[str, List[str]]]:
        """
        Fetch PubChem xrefs in batch mode by CID using PUG REST.

        Why this version:
          - Large batches can trigger PubChem 504s.
          - We *adaptively split* a failing batch into smaller chunks until it succeeds
            (down to `min_batch_size`) to keep overall retrieval time low and avoid
            losing xrefs entirely.

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

        # NOTE: PubChem xrefs endpoint does NOT reliably expose MeSH IDs.
        #       MeSH enrichment is handled via PUG-View (and optional MeSH RDF resolution).
        xref_types = "RegistryID,PubMedID,MIMID,GeneID,PatentID,DBURL,TaxonomyID"
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"

        batch_size = max(int(batch_size), 1)
        min_batch_size = max(int(min_batch_size), 1)

        logger.info(
            "Fetching PubChem xrefs for %d unique CIDs (batch_size=%d, min_batch_size=%d)",
            len(uniq_cids),
            batch_size,
            min_batch_size,
        )

        # Work queue of CID chunks; failing chunks will be split and re-queued.
        queue = deque([uniq_cids[i : i + batch_size] for i in range(0, len(uniq_cids), batch_size)])
        processed = 0

        while queue:
            batch = queue.popleft()
            processed += len(batch)

            payload = {"cid": ",".join(str(c) for c in batch)}
            success = False
            data: Dict[str, Any] = {}

            for attempt in range(1, int(max_retries) + 1):
                logger.info(
                    "PubChem xrefs batch (size=%d), processed ~%d/%d, attempt %d/%d",
                    len(batch),
                    min(processed, len(uniq_cids)),
                    len(uniq_cids),
                    attempt,
                    max_retries,
                )

                try:
                    resp = self.session.post(
                        f"{base_url}/xrefs/{xref_types}/JSON",
                        data=payload,
                        timeout=self.http_timeout,
                    )

                    # 404: no xrefs – treat as empty and stop retrying this chunk.
                    if resp.status_code == 404:
                        logger.warning("PubChem xrefs chunk (size=%d) returned 404 (no data).", len(batch))
                        success = True
                        data = {}
                        break

                    # 504: typically too-large batch; retry then split (below) if it keeps failing.
                    if resp.status_code == 504:
                        logger.warning("PubChem xrefs chunk (size=%d) got 504 (timeout).", len(batch))
                        if attempt < max_retries:
                            time.sleep(float(retry_sleep))
                        continue

                    resp.raise_for_status()
                    data = resp.json() if resp.content else {}
                    success = True
                    break

                except Exception as e:
                    logger.warning(
                        "PubChem xrefs chunk (size=%d) failed on attempt %d/%d: %s",
                        len(batch),
                        attempt,
                        max_retries,
                        e,
                    )
                    if attempt < max_retries:
                        time.sleep(float(retry_sleep))

            if not success:
                # Split & retry if chunk is still sizeable.
                if len(batch) > min_batch_size:
                    mid = len(batch) // 2
                    left = batch[:mid]
                    right = batch[mid:]
                    logger.warning(
                        "Splitting failing PubChem xrefs chunk of %d CIDs into %d + %d",
                        len(batch),
                        len(left),
                        len(right),
                    )
                    # Re-queue (left first so results appear earlier)
                    if right:
                        queue.appendleft(right)
                    if left:
                        queue.appendleft(left)
                else:
                    logger.error(
                        "PubChem xrefs chunk permanently failed (size=%d) after %d attempts; skipping these CIDs.",
                        len(batch),
                        max_retries,
                    )
                continue

            info_list = (
                data.get("InformationList", {}).get("Information", [])
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

            time.sleep(getattr(self, "pubchem_polite_sleep", 0.2))  # be polite to PubChem

        return results


    def enrich_compounds(
        self,
        input_name: str = "Compound_Properties_Processed.csv",
        output_name: str = "Compound_Properties_WithOntologies.csv",
        inchikey_col: str = "InChIKey",
        cid_col: str = "CompoundID",
        include_mesh_terms: bool = True,
        include_chemical_classes: bool = True,
        pugview_max_workers: int = 6,
        chebi_max_workers: int = 6,
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
                batch_size=getattr(self, "pubchem_xrefs_batch_size", 50),
                max_retries=getattr(self, "pubchem_xrefs_max_retries", 3),
                retry_sleep=getattr(self, "pubchem_xrefs_retry_sleep", 5.0),
                min_batch_size=getattr(self, "pubchem_xrefs_min_batch_size", 10),
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
            "MeSH_IDs": "MeSH_IDs",
            "DTXSID_IDs": "DTXSID_IDs",
            "KEGG_IDs": "KEGG_IDs",
            "HMDB_IDs": "HMDB_IDs",
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
            "MeSH_IDs": "MeSH_Sources",
            "DTXSID_IDs": "DTXSID_Sources",
            "KEGG_IDs": "KEGG_Sources",
            "HMDB_IDs": "HMDB_Sources",
        }

        # Track per-ontology, per-source counts for the run summary
        ontology_source_counts: Dict[str, Dict[str, int]] = {
            out_col: {"MyChem": 0, "PubChem": 0, "PUGView": 0}
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
        # 3b) Optional: MeSH terms + chemical classes (PubChem PUG-View + NLM MeSH + ChEBI parents)
        # ------------------------------------------------------
        cid_to_record: Dict[int, Dict[str, Any]] = {}
        cid_mesh_entry_terms: Dict[int, List[str]] = {}
        cid_mesh_pharm_class: Dict[int, List[str]] = {}
        cid_mesh_ref_ids: Dict[int, List[str]] = {}
        cid_mesh_embedded_ids: Dict[int, List[str]] = {}
        cid_chebi_embedded_ids: Dict[int, List[str]] = {}
        cid_pubchem_chem_class: Dict[int, Dict[str, List[str]]] = {}
        cid_pubchem_ontology_summary: Dict[int, str] = {}
        chebi_to_parents: Dict[str, Dict[str, List[str]]] = {}

        if (include_mesh_terms or include_chemical_classes) and cid_series is not None:
            uniq_cids = sorted({int(c) for c in cid_series.dropna().astype(int).tolist()})
            if uniq_cids:
                logger.info(
                    "Fetching PubChem PUG-View records for %d unique CIDs (max_workers=%d)…",
                    len(uniq_cids),
                    pugview_max_workers,
                )

                def _fetch_pug(cid: int):
                    return cid, self._pubchem_pug_view_record(cid)

                with ThreadPoolExecutor(max_workers=max(1, int(pugview_max_workers))) as ex:
                    futures = {ex.submit(_fetch_pug, cid): cid for cid in uniq_cids}
                    for fut in as_completed(futures):
                        cid, rec = fut.result()
                        if isinstance(rec, dict):
                            cid_to_record[cid] = rec

                # Parse PUG-View records (local CPU work)
                for cid, rec in cid_to_record.items():
                    try:
                        embedded = self._extract_ids_from_record_text(rec)
                        cid_mesh_embedded_ids[cid] = embedded.get("mesh_unique_ids", []) or []
                        cid_chebi_embedded_ids[cid] = embedded.get("chebi_ids", []) or []
                    except Exception:
                        cid_mesh_embedded_ids[cid] = []
                        cid_chebi_embedded_ids[cid] = []

                    if include_mesh_terms:
                        cid_mesh_entry_terms[cid] = self._extract_terms_by_toc_heading(rec, "MeSH Entry Terms")
                        cid_mesh_pharm_class[cid] = self._extract_terms_by_toc_heading(rec, "MeSH Pharmacological Classification")
                        refs = self._extract_pubchem_mesh_references(rec)
                        cid_mesh_ref_ids[cid] = [r.get("mesh_source_id", "") for r in refs if r.get("mesh_source_id")]

                    if include_chemical_classes:
                        class_secs = self._find_sections_heading_contains(rec, ("chemical classification", "chemical taxonomy"))
                        cid_pubchem_chem_class[cid] = self._extract_kv_from_sections(class_secs) if class_secs else {}
                        ont = self._extract_terms_by_toc_heading(rec, "Ontology Summary")
                        cid_pubchem_ontology_summary[cid] = (ont[0] if ont else "")

                # ChEBI parents (cached + parallel). We use both MyChem ChEBI IDs and any ChEBI tokens embedded in PubChem.
                if include_chemical_classes:
                    chebi_ids: Set[str] = set()

                    if "ChEBI_IDs" in df.columns:
                        for s in df["ChEBI_IDs"].fillna("").astype(str).tolist():
                            for tok in s.split("|"):
                                tok = tok.strip()
                                if tok:
                                    chebi_ids.add(tok)

                    for lst in cid_chebi_embedded_ids.values():
                        for tok in lst or []:
                            tok = str(tok).strip()
                            if tok:
                                chebi_ids.add(tok)

                    chebi_ids_list = sorted(chebi_ids)
                    if chebi_ids_list:
                        logger.info(
                            "Fetching ChEBI parent classes for %d unique ChEBI IDs (max_workers=%d)…",
                            len(chebi_ids_list),
                            chebi_max_workers,
                        )

                        def _fetch_chebi(cid: str):
                            parents = self._chebi_get_parent_nodes(cid)
                            p_ids = sorted({p.get("chebiId") for p in parents if p.get("chebiId")})
                            p_names = sorted({p.get("chebiName") for p in parents if p.get("chebiName")})
                            return cid, [str(x) for x in p_ids if x], [str(x) for x in p_names if x]

                        with ThreadPoolExecutor(max_workers=max(1, int(chebi_max_workers))) as ex:
                            futures = {ex.submit(_fetch_chebi, c): c for c in chebi_ids_list}
                            for fut in as_completed(futures):
                                c, p_ids, p_names = fut.result()
                                chebi_to_parents[c] = {"ids": p_ids, "names": p_names}

        # Create new columns (row-wise) without extra HTTP calls
        if include_mesh_terms or include_chemical_classes:
            mesh_entry_out: List[str] = []
            mesh_pharm_out: List[str] = []
            mesh_embedded_ids_out: List[str] = []
            mesh_all_ids_out: List[str] = []
            mesh_terms_sources_out: List[str] = []

            pubchem_ontology_summary_out: List[str] = []
            pubchem_chemclass_json_out: List[str] = []
            pubchem_kingdom_out: List[str] = []
            pubchem_superclass_out: List[str] = []
            pubchem_class_out: List[str] = []
            pubchem_subclass_out: List[str] = []
            pubchem_direct_parent_out: List[str] = []
            pubchem_alt_parent_out: List[str] = []
            pubchem_framework_out: List[str] = []

            chebi_parent_ids_out: List[str] = []
            chebi_parent_names_out: List[str] = []
            chemical_classes_sources_out: List[str] = []

            for _, row in df.iterrows():
                cid_val: Optional[int] = None
                if cid_series is not None:
                    cid_raw = row.get(cid_col)
                    if pd.notna(cid_raw):
                        try:
                            cid_val = int(cid_raw)
                        except Exception:
                            cid_val = None

                # Existing MeSH IDs (from MyChem/PubChem xrefs) – keep as part of completeness union
                existing_mesh_ids: Set[str] = set()
                existing_mesh_raw = row.get("MeSH_IDs", "")
                if pd.notna(existing_mesh_raw):
                    for tok in str(existing_mesh_raw).split("|"):
                        tok = tok.strip()
                        if tok:
                            existing_mesh_ids.add(tok)

                entry_terms = cid_mesh_entry_terms.get(cid_val, []) if cid_val is not None else []
                pharm_terms = cid_mesh_pharm_class.get(cid_val, []) if cid_val is not None else []
                embedded_mesh_ids = cid_mesh_embedded_ids.get(cid_val, []) if cid_val is not None else []
                pubchem_ref_ids = cid_mesh_ref_ids.get(cid_val, []) if cid_val is not None else []
                # Union all MeSH IDs we know about (completeness)
                all_mesh_ids = set(existing_mesh_ids)
                all_mesh_ids.update([x for x in embedded_mesh_ids if x])
                all_mesh_ids.update([x for x in pubchem_ref_ids if x])

                # Sources
                mesh_sources: Set[str] = set()
                if existing_mesh_ids:
                    mesh_sources.add("MyChem/PubChemXrefs")
                if (entry_terms or pharm_terms or embedded_mesh_ids or pubchem_ref_ids) and cid_val is not None:
                    mesh_sources.add("PubChemPUGView")

                mesh_entry_out.append("|".join(entry_terms) if include_mesh_terms and entry_terms else "")
                mesh_pharm_out.append("|".join(pharm_terms) if include_mesh_terms and pharm_terms else "")
                mesh_embedded_ids_out.append("|".join(sorted(set(embedded_mesh_ids))) if include_mesh_terms and embedded_mesh_ids else "")
                mesh_all_ids_out.append("|".join(sorted(all_mesh_ids)) if include_mesh_terms and all_mesh_ids else "")
                mesh_terms_sources_out.append(";".join(sorted(mesh_sources)) if include_mesh_terms and mesh_sources else "")

                # Chemical classes
                chem_sources: Set[str] = set()
                chem_dict: Dict[str, List[str]] = {}
                ont_sum = ""
                if cid_val is not None and include_chemical_classes:
                    chem_dict = cid_pubchem_chem_class.get(cid_val, {}) or {}
                    ont_sum = cid_pubchem_ontology_summary.get(cid_val, "") or ""
                    if chem_dict or ont_sum:
                        chem_sources.add("PubChemPUGView")

                # Extract commonly used ClassyFire-like keys if present
                def _pick(keys: Tuple[str, ...]) -> str:
                    for k in keys:
                        if k in chem_dict and chem_dict.get(k):
                            return "|".join([str(x) for x in chem_dict.get(k) if x])
                    return ""

                pubchem_ontology_summary_out.append(ont_sum if include_chemical_classes else "")
                pubchem_chemclass_json_out.append(json.dumps(chem_dict, ensure_ascii=False) if include_chemical_classes and chem_dict else "")
                pubchem_kingdom_out.append(_pick(("Kingdom", "kingdom")))
                pubchem_superclass_out.append(_pick(("Superclass", "superclass")))
                pubchem_class_out.append(_pick(("Class", "class")))
                pubchem_subclass_out.append(_pick(("Subclass", "subclass")))
                pubchem_direct_parent_out.append(_pick(("Direct Parent", "Direct parent", "DirectParent")))
                pubchem_alt_parent_out.append(_pick(("Alternative Parent", "Alternative parent", "AlternativeParent")))
                pubchem_framework_out.append(_pick(("Molecular Framework", "Molecular framework", "MolecularFramework")))

                # ChEBI parents aggregated
                parent_ids: Set[str] = set()
                parent_names: Set[str] = set()

                chebi_ids_for_row: Set[str] = set()
                chebi_raw = row.get("ChEBI_IDs", "")
                if pd.notna(chebi_raw):
                    for tok in str(chebi_raw).split("|"):
                        tok = tok.strip()
                        if tok:
                            chebi_ids_for_row.add(tok)

                if cid_val is not None:
                    for tok in cid_chebi_embedded_ids.get(cid_val, []) or []:
                        tok = str(tok).strip()
                        if tok:
                            chebi_ids_for_row.add(tok)

                for chebi_id in chebi_ids_for_row:
                    p = chebi_to_parents.get(chebi_id)
                    if p:
                        parent_ids.update(p.get("ids", []) or [])
                        parent_names.update(p.get("names", []) or [])

                if chebi_ids_for_row:
                    chem_sources.add("ChEBI")

                chebi_parent_ids_out.append("|".join(sorted(parent_ids)) if include_chemical_classes and parent_ids else "")
                chebi_parent_names_out.append("|".join(sorted(parent_names)) if include_chemical_classes and parent_names else "")
                chemical_classes_sources_out.append(";".join(sorted(chem_sources)) if include_chemical_classes and chem_sources else "")

            # Attach new columns
            if include_mesh_terms:
                df["MeSH_Entry_Terms"] = mesh_entry_out
                df["MeSH_Pharmacological_Classification"] = mesh_pharm_out
                df["MeSH_Embedded_IDs"] = mesh_embedded_ids_out
                df["MeSH_All_IDs"] = mesh_all_ids_out
                df["MeSH_Terms_Sources"] = mesh_terms_sources_out

            # Backfill legacy MeSH_IDs + MeSH_Sources to include PUG-View / MeSH RDF resolution.
            # This makes downstream graph linking easier while keeping provenance explicit.
            if "MeSH_IDs" in df.columns and "MeSH_All_IDs" in df.columns and "MeSH_Sources" in df.columns:
                merged_ids: List[str] = []
                merged_sources: List[str] = []
                pugview_rows = 0

                for existing_ids, all_ids, srcs in zip(
                    df["MeSH_IDs"].fillna(""),
                    df["MeSH_All_IDs"].fillna(""),
                    df["MeSH_Sources"].fillna(""),
                ):
                    ids_set = set([x for x in str(existing_ids).split("|") if x]) | set([x for x in str(all_ids).split("|") if x])
                    merged_ids.append("|".join(sorted(ids_set)) if ids_set else "")

                    src_set = set([x for x in str(srcs).split(";") if x])
                    if str(all_ids).strip():
                        src_set.add("PUGView")
                        pugview_rows += 1
                    merged_sources.append(";".join(sorted(src_set)) if src_set else "")

                df["MeSH_IDs"] = merged_ids
                df["MeSH_Sources"] = merged_sources

                # Update per-source counters for MeSH only (other ontologies remain MyChem/PubChem only)
                ontology_source_counts["MeSH_IDs"]["PUGView"] = int(pugview_rows)

            if include_chemical_classes:
                df["PubChem_Ontology_Summary"] = pubchem_ontology_summary_out
                df["PubChem_Chemical_Classification_JSON"] = pubchem_chemclass_json_out
                df["PubChem_ChemicalClass_Kingdom"] = pubchem_kingdom_out
                df["PubChem_ChemicalClass_Superclass"] = pubchem_superclass_out
                df["PubChem_ChemicalClass_Class"] = pubchem_class_out
                df["PubChem_ChemicalClass_Subclass"] = pubchem_subclass_out
                df["PubChem_ChemicalClass_DirectParent"] = pubchem_direct_parent_out
                df["PubChem_ChemicalClass_AlternativeParent"] = pubchem_alt_parent_out
                df["PubChem_ChemicalClass_MolecularFramework"] = pubchem_framework_out
                df["ChEBI_Parent_IDs"] = chebi_parent_ids_out
                df["ChEBI_Parent_Names"] = chebi_parent_names_out
                df["Chemical_Classes_Sources"] = chemical_classes_sources_out

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
            fields_to_check = list(output_fields.keys())
            if include_mesh_terms:
                fields_to_check += [
                    "MeSH_Entry_Terms",
                    "MeSH_Pharmacological_Classification",
                    "MeSH_All_IDs",
                ]
            if include_chemical_classes:
                fields_to_check += [
                    "PubChem_Chemical_Classification_JSON",
                    "PubChem_Ontology_Summary",
                    "ChEBI_Parent_IDs",
                ]

            has_any = any(bool(str(row.get(col, "")).strip()) for col in fields_to_check)
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
                "  %s: %d  (%.2f%%) [MyChem=%d, PubChem=%d, PUGView=%d]",
                k,
                v,
                ontology_percent[k],
                ontology_source_counts[k]["MyChem"],
                ontology_source_counts[k]["PubChem"],
                ontology_source_counts[k].get("PUGView", 0)
            )

        # Build final summary dict
        retrieved_fields = list(output_fields.keys())
        if include_mesh_terms:
            retrieved_fields += [
                "MeSH_Entry_Terms",
                "MeSH_Pharmacological_Classification",
                "MeSH_Embedded_IDs",
                "MeSH_All_IDs",
                "MeSH_Terms_Sources",
            ]
        if include_chemical_classes:
            retrieved_fields += [
                "PubChem_Ontology_Summary",
                "PubChem_Chemical_Classification_JSON",
                "PubChem_ChemicalClass_Kingdom",
                "PubChem_ChemicalClass_Superclass",
                "PubChem_ChemicalClass_Class",
                "PubChem_ChemicalClass_Subclass",
                "PubChem_ChemicalClass_DirectParent",
                "PubChem_ChemicalClass_AlternativeParent",
                "PubChem_ChemicalClass_MolecularFramework",
                "ChEBI_Parent_IDs",
                "ChEBI_Parent_Names",
                "Chemical_Classes_Sources",
            ]

        run_summary = {
            "total_compounds_processed": int(num_total),
            "total_nodes_processed": int(num_total),  # generic key for tooling/tests
            "total_unique_inchikeys": int(num_total),
            "total_successful_hits": int(num_total - num_missing),
            "total_successful_hits_any": int(num_total - num_missing),  # generic key for tooling/tests
            "total_missing_entries": int(num_missing),
            "missing_inchikeys": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "ontology_source_counts": ontology_source_counts,  # <-- NEW
            "retrieved_fields": retrieved_fields,
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
            "batch_size_pubchem": int(getattr(self, "pubchem_xrefs_batch_size", 50)),
            "source": "MyChem.info(BioThings);PubChem PUG REST",
        }

        summary_path = f"{self.data_dir}/Compound_Ontology_RunSummary.json"
        self._save_run_summary(run_summary, summary_path)

        # ------------------------------------------------------
        # 5) Save enriched version
        # ------------------------------------------------------
        df["OntologyEnrichedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.MyChem+PubChem+PUGView+MeSH+ChEBI/1.2"
        df["OntologySources"] = "MyChem.info;PubChem;PubChemPUGView;NLMMeSH;ChEBI"

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
                "pathway.kegg",
                "taxid",
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
                # KEGG pathways (IDs)
                # -------------------------------
                kegg_ids: List[str] = []
                kegg = None
                if isinstance(pathway, dict):
                    kegg = pathway.get("kegg")

                if isinstance(kegg, dict):
                    kegg = [kegg]

                if isinstance(kegg, list):
                    for pw in kegg:
                        if not isinstance(pw, dict):
                            continue
                        pid = pw.get("id") or pw.get("pathway_id")
                        if pid:
                            kegg_ids.append(str(pid))

                if kegg_ids:
                    out["KEGG"] = self._uniq(kegg_ids)

                # -------------------------------
                # Taxonomy (NCBITaxon)
                # -------------------------------
                taxid = doc.get("taxid") or doc.get("taxon") or doc.get("taxon_id")
                if isinstance(taxid, (str, int)):
                    out["TaxonomyID"] = [str(taxid)]
                elif isinstance(taxid, list):
                    out["TaxonomyID"] = self._uniq([str(x) for x in taxid if x])

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
            "KEGG_Pathway_IDs": "KEGG",
            "NCBITaxon_ID": "TaxonomyID",
            "OMIM_IDs": "OMIM",
            # Reserve disease/phenotype ontologies
            # "DOID_IDs": "DOID",    # not populated (requires external KPs)
            # "MONDO_IDs": "MONDO",  # not populated
            # "HPO_Terms": "HPO",    # not populated
        }

        retrieved_fields = list(output_fields.keys())

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
                    "KEGG_Pathway_IDs",
                    "NCBITaxon_ID",
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
            "KEGG_Pathway_IDs": "MyGene.info (pathway.kegg)",
            "NCBITaxon_ID": "MyGene.info (taxid)",
            "OMIM_IDs": "MyGene.info (omim)",
            # "DOID_IDs": "Not populated – requires disease-centric KPs (e.g. Monarch/MyDisease)",
            # "MONDO_IDs": "Not populated – requires disease-centric KPs (e.g. ClinGen/Monarch)",
            # "HPO_Terms": "Not populated – requires phenotype-centric KPs (e.g. HPO/Monarch)",
        }

        # Build and save JSON summary
        run_summary = {
            "total_genes_processed": int(num_total),
            "total_nodes_processed": int(num_total),  # generic key for tooling/tests
            "total_unique_gene_ids": int(num_total),
            "total_successful_hits": int(num_total - num_missing),
            "total_successful_hits_any": int(num_total - num_missing),  # generic key for tooling/tests
            "total_missing_entries": int(num_missing),
            "missing_gene_ids": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "retrieved_fields": retrieved_fields,
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
                "pathway.kegg",
                "taxid",
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
        fields = "accession,ec,go_p,go_f,go_c," + ",".join(xref_field_list)

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
                while len(cols) < 5 + len(xref_field_list):
                    cols.append("")

                acc = cols[0].strip()
                if not acc:
                    continue

                ec_raw = cols[1]
                bp_raw = cols[2]
                mf_raw = cols[3]
                cc_raw = cols[4]

                def _split_ec_field(raw: str) -> List[str]:
                    """Extract EC numbers from UniProt EC field (TSV)."""
                    if not raw:
                        return []
                    import re
                    toks: List[str] = []
                    for piece in str(raw).replace(",", ";").split(";"):
                        piece = piece.strip()
                        if not piece:
                            continue
                        m = re.search(r"(\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-))", piece)
                        if m:
                            toks.append(m.group(1))
                    return toks

                ann: Dict[str, List[str]] = {
                    "EC": self._uniq(_split_ec_field(ec_raw)),
                    "GO_BP": self._uniq(_split_go_field(bp_raw)),
                    "GO_MF": self._uniq(_split_go_field(mf_raw)),
                    "GO_CC": self._uniq(_split_go_field(cc_raw)),
                }

                # Parse all configured xref_* columns into labelled lists
                for offset, field_name in enumerate(xref_field_list):
                    raw_val = cols[5 + offset]
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
        ec_numbers: List[str] = []

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
            ec_numbers.append(_join(ann.get("EC", [])))

            for label, col_name in xref_output_map.items():
                xref_values_lists[col_name].append(_join(ann.get(label, [])))

        df["GO_BP_Protein_IDs"] = go_bp_ids
        df["GO_MF_Protein_IDs"] = go_mf_ids
        df["GO_CC_Protein_IDs"] = go_cc_ids
        df["EC_Numbers"] = ec_numbers

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
            "EC_Numbers",
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
            "EC_Numbers",
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
            "total_nodes_processed": int(num_rows),  # generic key for tooling/tests
            "total_unique_input_accessions": int(len(raw_ids)),
            "total_successful_hits": int(num_success),
            "total_successful_hits_any": int(num_success),  # generic key for tooling/tests
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
        endpoint_cols: Optional[List[str]] = None,
        ols_max_workers: int = 4,
    ) -> Optional[pd.DataFrame]:
        """
        Enrich BioAssay nodes with assay ontologies via EBI OLS.

        Covers the schema requirements from FigureA/FigureB:
          - Assay ontologies: BAO + OBI + EFO
          - Endpoint terms: BAO/EFO (derived from endpoint-like columns)
          - AssayFormat (biochemical/cellular/microsomes) via lightweight rules

        Output columns added (per ontology = BAO, OBI, EFO):
            *_ID_String, *_IRI_String, *_Label_String, *_StringScore
            *_ID_Keyword, *_IRI_Keyword, *_Label_Keyword, *_KeywordScore
            *_ID, *_IRI, *_Label, *_Score   (combined / main mapping)

        Additional columns:
            EndpointTerm_Raw
            EndpointTerm_BAO_IDs, EndpointTerm_BAO_Labels
            EndpointTerm_EFO_IDs, EndpointTerm_EFO_Labels
            EndpointTerm_Sources
            AssayFormat, AssayFormat_Confidence, AssayFormat_Method

        A JSON summary is written to:
            Data/Nodes/Assay_Ontology_RunSummary.json
        """
        input_path = f"{self.data_dir}/{input_name}"

        try:
            df = pd.read_csv(input_path, low_memory=False)
        except FileNotFoundError:
            logger.info("Assay input file not found: %s (skipping)", input_path)
            return None

        if name_cols is None:
            name_cols = ["AssayName", "AssayActivityName", "AssayType", "AssayDescription"]

        if endpoint_cols is None:
            endpoint_cols = [
                "AssayActivityName",
                "AssayEndpoint",
                "EndpointName",
                "ReadoutType",
                "AssayResultName",
                "AssayDescription",
            ]

        # -----------------------------
        # Helpers
        # -----------------------------
        def _iri_to_curie(iri: str) -> str:
            if not iri:
                return ""
            local = iri.rsplit("/", 1)[-1]
            if "#" in local:
                local = local.split("#")[-1]
            if "_" in local:
                prefix, rest = local.split("_", 1)
                return f"{prefix}:{rest}"
            return local

        def _infer_assay_format(text: str) -> Tuple[str, float]:
            """
            Very lightweight normalizer for AssayFormat node.
            Returns (format_label, confidence).
            """
            t = (text or "").lower()
            if not t:
                return ("", 0.0)

            # microsomes / S9 / liver fractions
            if any(k in t for k in ["microsome", "microsomes", "s9 fraction", "s9 ", "liver microsome", "hepatic microsome"]):
                return ("microsomes", 0.9)

            # cell/tissue-based
            if any(k in t for k in ["cell line", "cell-based", "cell based", "cellular", "hepg2", "primary hepatocyte", "in vivo", "tissue"]):
                return ("cellular", 0.8)

            # biochemical / enzyme / binding
            if any(k in t for k in ["enzyme", "activity assay", "inhibition", "binding", "biochemical", "kinase", "assay buffer", "substrate"]):
                return ("biochemical", 0.75)

            return ("", 0.0)

        def _split_endpoints(raw: str) -> List[str]:
            if not raw:
                return []
            # split common delimiters; keep short phrases
            parts = re.split(r"[;|,/]+", str(raw))
            out = []
            for p in parts:
                p = p.strip()
                if not p:
                    continue
                # drop extremely generic tokens
                if p.lower() in {"other", "na", "n/a", "none"}:
                    continue
                out.append(p[:120])
            return out

        # -----------------------------
        # Build search texts (assay-level) + endpoint tokens
        # -----------------------------
        # Use existing lightweight builder (prioritizes informative columns)
        search_texts: List[str] = [self._build_search_text(row) for _, row in df.iterrows()]

        endpoint_tokens_per_row: List[List[str]] = []
        for _, row in df.iterrows():
            toks: List[str] = []
            for c in endpoint_cols:
                if c in row and pd.notna(row[c]):
                    toks.extend(_split_endpoints(str(row[c])))
            # de-dup per row, keep order
            seen = set()
            uniq = []
            for t in toks:
                tl = t.lower()
                if tl in seen:
                    continue
                seen.add(tl)
                uniq.append(t)
            endpoint_tokens_per_row.append(uniq)

        unique_search_texts = sorted({t for t in search_texts if t})
        unique_endpoints = sorted({t for toks in endpoint_tokens_per_row for t in toks if t})

        logger.info(
            "Enriching %d assays (unique assay texts=%d, unique endpoints=%d) using OLS BAO/OBI/EFO",
            len(df),
            len(unique_search_texts),
            len(unique_endpoints),
        )

        # -----------------------------
        # OLS lookup (multi-score)
        # -----------------------------
        def _ols_lookup_multi(term_clean: str, ontology: str) -> Dict[str, Tuple[str, str, float]]:
            """
            Returns:
              {"string": (iri,label,score), "keyword": (...), "combined": (...)}
            """
            if not term_clean:
                return {"string": ("", "", 0.0), "keyword": ("", "", 0.0), "combined": ("", "", 0.0)}

            url = f"{self.ols_base_url}/search"
            params = {
                "q": term_clean,
                "ontology": ontology,
                "rows": 10,
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
                    if attempt == 3:
                        logger.debug("OLS failed for '%s' (%s): %s", term_clean, ontology, e)
                        return {"string": ("", "", 0.0), "keyword": ("", "", 0.0), "combined": ("", "", 0.0)}
                    time.sleep(0.4 * attempt)

            docs = (data or {}).get("response", {}).get("docs", [])
            if not docs:
                return {"string": ("", "", 0.0), "keyword": ("", "", 0.0), "combined": ("", "", 0.0)}

            term_lower = term_clean.lower()
            src_keywords = self._extract_keywords(term_clean)

            best_string = ("", "", 0.0)
            best_keyword = ("", "", 0.0)
            best_combined = ("", "", 0.0)

            for d in docs:
                label = (d.get("label") or term_clean).strip()
                label_lower = label.lower()

                synonyms = d.get("synonym") or []
                if isinstance(synonyms, str):
                    synonyms = [synonyms]

                short_form = d.get("short_form") or []
                if isinstance(short_form, str):
                    short_form = [short_form]

                # 1) String similarity
                if label_lower == term_lower or any(str(s).strip().lower() == term_lower for s in synonyms):
                    string_score = 1.0
                else:
                    string_score = SequenceMatcher(None, term_lower, label_lower).ratio()

                # 2) Keyword overlap
                if src_keywords:
                    target_text = " ".join([label] + [str(s) for s in synonyms] + [str(sf) for sf in short_form])
                    tgt_keywords = self._extract_keywords(target_text)
                    keyword_score = (len(src_keywords & tgt_keywords) / len(src_keywords)) if tgt_keywords else 0.0
                else:
                    keyword_score = 0.0

                # 3) Combined
                alpha = 0.7
                combined_score = alpha * string_score + (1.0 - alpha) * keyword_score

                iri = d.get("iri") or ""

                if string_score > best_string[2]:
                    best_string = (iri, label, float(string_score))
                if keyword_score > best_keyword[2]:
                    best_keyword = (iri, label, float(keyword_score))
                if combined_score > best_combined[2]:
                    best_combined = (iri, label, float(combined_score))

            # mild floor to avoid very weak mappings
            min_combined = 0.35 if ontology in {"bao", "obi", "efo"} else 0.30
            if best_combined[2] < min_combined:
                best_combined = ("", "", 0.0)

            if self.OLS_SLEEP_SECONDS and self.OLS_SLEEP_SECONDS > 0:
                time.sleep(self.OLS_SLEEP_SECONDS)

            return {"string": best_string, "keyword": best_keyword, "combined": best_combined}

        def _bulk_lookup(terms: List[str], ontology: str) -> Dict[str, Dict[str, Tuple[str, str, float]]]:
            if not terms:
                return {}
            out: Dict[str, Dict[str, Tuple[str, str, float]]] = {}
            if ols_max_workers <= 1:
                for t in terms:
                    out[t] = _ols_lookup_multi(t, ontology)
                return out

            from concurrent.futures import ThreadPoolExecutor, as_completed

            with ThreadPoolExecutor(max_workers=ols_max_workers) as ex:
                futs = {ex.submit(_ols_lookup_multi, t, ontology): t for t in terms}
                for fut in as_completed(futs):
                    t = futs[fut]
                    try:
                        out[t] = fut.result()
                    except Exception:  # pragma: no cover
                        out[t] = {"string": ("", "", 0.0), "keyword": ("", "", 0.0), "combined": ("", "", 0.0)}
            return out

        # -----------------------------
        # 1) Assay-level lookups
        # -----------------------------
        assay_maps = {
            "bao": _bulk_lookup(unique_search_texts, "bao"),
            "obi": _bulk_lookup(unique_search_texts, "obi"),
            "efo": _bulk_lookup(unique_search_texts, "efo"),
        }

        # -----------------------------
        # 2) Endpoint term lookups (combined only)
        # -----------------------------
        endpoint_maps = {
            "bao": _bulk_lookup(unique_endpoints, "bao"),
            "efo": _bulk_lookup(unique_endpoints, "efo"),
        }

        # -----------------------------
        # Attach columns (per ontology)
        # -----------------------------
        def _attach_ontology_cols(prefix: str, ontology_key: str) -> None:
            id_s, iri_s, lab_s, sc_s = [], [], [], []
            id_k, iri_k, lab_k, sc_k = [], [], [], []
            ids, iris, labs, sc = [], [], [], []

            mapping = assay_maps[ontology_key]

            for text in search_texts:
                res = mapping.get(text) if text else None
                if not res:
                    res = {"string": ("", "", 0.0), "keyword": ("", "", 0.0), "combined": ("", "", 0.0)}

                i_s, l_s, s_s = res["string"]
                i_k, l_k, s_k = res["keyword"]
                i_c, l_c, s_c = res["combined"]

                iri_s.append(i_s); lab_s.append(l_s); sc_s.append(float(s_s)); id_s.append(_iri_to_curie(i_s))
                iri_k.append(i_k); lab_k.append(l_k); sc_k.append(float(s_k)); id_k.append(_iri_to_curie(i_k))
                iris.append(i_c); labs.append(l_c); sc.append(float(s_c)); ids.append(_iri_to_curie(i_c))

            df[f"{prefix}_ID_String"] = id_s
            df[f"{prefix}_IRI_String"] = iri_s
            df[f"{prefix}_Label_String"] = lab_s
            df[f"{prefix}_StringScore"] = sc_s

            df[f"{prefix}_ID_Keyword"] = id_k
            df[f"{prefix}_IRI_Keyword"] = iri_k
            df[f"{prefix}_Label_Keyword"] = lab_k
            df[f"{prefix}_KeywordScore"] = sc_k

            df[f"{prefix}_ID"] = ids
            df[f"{prefix}_IRI"] = iris
            df[f"{prefix}_Label"] = labs
            df[f"{prefix}_Score"] = sc

        _attach_ontology_cols("BAO", "bao")
        _attach_ontology_cols("OBI", "obi")
        _attach_ontology_cols("EFO", "efo")

        # -----------------------------
        # Endpoint term columns (combined only)
        # -----------------------------
        endpoint_raw_col: List[str] = []
        endpoint_bao_ids: List[str] = []
        endpoint_bao_labels: List[str] = []
        endpoint_efo_ids: List[str] = []
        endpoint_efo_labels: List[str] = []

        def _join_unique(xs: List[str]) -> str:
            xs = [x for x in xs if x]
            return "|".join(sorted(set(xs))) if xs else ""

        for toks in endpoint_tokens_per_row:
            endpoint_raw_col.append("|".join(toks) if toks else "")

            bao_ids, bao_labs = [], []
            efo_ids, efo_labs = [], []

            for t in toks:
                r_bao = endpoint_maps["bao"].get(t, {}).get("combined", ("", "", 0.0))
                r_efo = endpoint_maps["efo"].get(t, {}).get("combined", ("", "", 0.0))

                iri_b, lab_b, _ = r_bao
                iri_e, lab_e, _ = r_efo

                cid_b = _iri_to_curie(iri_b)
                cid_e = _iri_to_curie(iri_e)

                if cid_b:
                    bao_ids.append(cid_b)
                    bao_labs.append(lab_b)
                if cid_e:
                    efo_ids.append(cid_e)
                    efo_labs.append(lab_e)

            endpoint_bao_ids.append(_join_unique(bao_ids))
            endpoint_bao_labels.append(_join_unique(bao_labs))
            endpoint_efo_ids.append(_join_unique(efo_ids))
            endpoint_efo_labels.append(_join_unique(efo_labs))

        df["EndpointTerm_Raw"] = endpoint_raw_col
        df["EndpointTerm_BAO_IDs"] = endpoint_bao_ids
        df["EndpointTerm_BAO_Labels"] = endpoint_bao_labels
        df["EndpointTerm_EFO_IDs"] = endpoint_efo_ids
        df["EndpointTerm_EFO_Labels"] = endpoint_efo_labels
        df["EndpointTerm_Sources"] = "OLS(BAO;EFO)"

        # -----------------------------
        # AssayFormat columns (rule-based)
        # -----------------------------
        fmt_vals: List[str] = []
        fmt_conf: List[float] = []
        fmt_method: List[str] = []

        for text in search_texts:
            fmt, conf = _infer_assay_format(text)
            fmt_vals.append(fmt)
            fmt_conf.append(conf)
            fmt_method.append("rule-based")

        df["AssayFormat"] = fmt_vals
        df["AssayFormat_Confidence"] = fmt_conf
        df["AssayFormat_Method"] = fmt_method

        # -----------------------------
        # Coverage summary (combined mappings)
        # -----------------------------
        ontology_cols = [
            "BAO_ID",
            "OBI_ID",
            "EFO_ID",
            "EndpointTerm_BAO_IDs",
            "EndpointTerm_EFO_IDs",
            "AssayFormat",
        ]

        ontology_counts: Dict[str, int] = {c: self._count_non_empty_column(df[c]) for c in ontology_cols}

        num_rows = len(df)
        ontology_percent = {
            k: (v / num_rows * 100.0 if num_rows > 0 else 0.0)
            for k, v in ontology_counts.items()
        }

        # Identify assays with no ontology at all (BAO/OBI/EFO empty AND no endpoint term)
        missing_all: List[int] = []
        id_col = "AssayID" if "AssayID" in df.columns else None

        for _, row in df.iterrows():
            has_any = any(pd.notna(row[c]) and str(row[c]).strip() != "" for c in ["BAO_ID", "OBI_ID", "EFO_ID", "EndpointTerm_BAO_IDs", "EndpointTerm_EFO_IDs"])
            if not has_any and id_col is not None:
                aid = row.get(id_col)
                if pd.notna(aid):
                    try:
                        missing_all.append(int(aid))
                    except Exception:
                        pass

        num_missing = len(missing_all)
        num_success = num_rows - num_missing

        logger.info("=== Assay Ontology Retrieval Summary ===")
        logger.info("Total assays processed: %d", num_rows)
        logger.info("Successful hits (any ontology/endpoint): %d", num_success)
        logger.info("Missing entries (no BAO/OBI/EFO/endpoint at all): %d", num_missing)
        for k, v in ontology_counts.items():
            logger.info("  %s: %d  (%.2f%%)", k, v, ontology_percent[k])

        run_summary = {
            "total_assays_processed": int(num_rows),
            "total_nodes_processed": int(num_rows),
            "total_successful_hits_any": int(num_success),
            "total_missing_entries": int(num_missing),
            "missing_assay_ids": missing_all,
            "ontology_counts": ontology_counts,
            "ontology_coverage_percent": ontology_percent,
            "retrieved_fields": (
                # BAO
                ["BAO_ID_String", "BAO_IRI_String", "BAO_Label_String", "BAO_StringScore",
                 "BAO_ID_Keyword", "BAO_IRI_Keyword", "BAO_Label_Keyword", "BAO_KeywordScore",
                 "BAO_ID", "BAO_IRI", "BAO_Label", "BAO_Score"]
                # OBI
                + ["OBI_ID_String", "OBI_IRI_String", "OBI_Label_String", "OBI_StringScore",
                   "OBI_ID_Keyword", "OBI_IRI_Keyword", "OBI_Label_Keyword", "OBI_KeywordScore",
                   "OBI_ID", "OBI_IRI", "OBI_Label", "OBI_Score"]
                # EFO
                + ["EFO_ID_String", "EFO_IRI_String", "EFO_Label_String", "EFO_StringScore",
                   "EFO_ID_Keyword", "EFO_IRI_Keyword", "EFO_Label_Keyword", "EFO_KeywordScore",
                   "EFO_ID", "EFO_IRI", "EFO_Label", "EFO_Score"]
                # endpoints + format
                + ["EndpointTerm_Raw", "EndpointTerm_BAO_IDs", "EndpointTerm_BAO_Labels",
                   "EndpointTerm_EFO_IDs", "EndpointTerm_EFO_Labels", "EndpointTerm_Sources",
                   "AssayFormat", "AssayFormat_Confidence", "AssayFormat_Method"]
            ),
            "input_file": input_path,
            "output_file": f"{self.data_dir}/{output_name}",
            "name_columns_used": name_cols,
            "endpoint_columns_used": endpoint_cols,
            "ols_ontologies": ["bao", "obi", "efo"],
            "ols_base_url": self.ols_base_url,
            "created_at": datetime.datetime.utcnow().isoformat() + "Z",
            "ols_sleep_seconds": self.OLS_SLEEP_SECONDS,
            "ols_max_workers": int(ols_max_workers),
            "source": "EBI OLS (BAO;OBI;EFO) + rule-based AssayFormat",
        }

        summary_path = f"{self.data_dir}/Assay_Ontology_RunSummary.json"
        self._save_run_summary(run_summary, summary_path)

        # Save enriched version
        enrichment_time = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnrichedAt"] = enrichment_time
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.AssaysOLS/2.0"
        df["OntologySources"] = "OLS(BAO;OBI;EFO);RuleBased(AssayFormat)"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved assay-enriched table to %s", out_path)

        return df

    # ------------------------------------------------------------------
    # Auxiliary schema nodes (FigureB): ExperimentalContext, Unit, Taxon, Tissue, CellLine, AssayEndpoint
    # These are optional: the methods auto-skip if the expected input file doesn't exist.
    # ------------------------------------------------------------------
    def _detect_first_column(self, df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
        for c in candidates:
            if c in df.columns:
                return c
        return None

    def _ols_best_match(
        self,
        term_clean: str,
        ontology: str,
        min_combined: float = 0.35,
        rows: int = 10,
    ) -> Tuple[str, str, str, float]:
        """
        Best-effort OLS matcher for controlled terms.

        Returns:
            (curie, iri, label, combined_score)
        """
        if not term_clean:
            return ("", "", "", 0.0)

        url = f"{self.ols_base_url}/search"
        params = {
            "q": term_clean,
            "ontology": ontology,
            "rows": rows,
            "queryFields": "label,synonym,short_form",
        }

        try:
            resp = self.session.get(url, params=params, timeout=self.http_timeout)
            resp.raise_for_status()
            data = resp.json()
        except Exception:  # pragma: no cover
            return ("", "", "", 0.0)

        docs = (data or {}).get("response", {}).get("docs", [])
        if not docs:
            return ("", "", "", 0.0)

        term_lower = term_clean.lower()
        src_keywords = self._extract_keywords(term_clean)

        best = ("", "", "", 0.0)  # curie, iri, label, score

        def iri_to_curie(iri: str) -> str:
            if not iri:
                return ""
            local = iri.rsplit("/", 1)[-1]
            if "#" in local:
                local = local.split("#")[-1]
            if "_" in local:
                prefix, rest = local.split("_", 1)
                return f"{prefix}:{rest}"
            return local

        for d in docs:
            label = (d.get("label") or term_clean).strip()
            label_lower = label.lower()

            synonyms = d.get("synonym") or []
            if isinstance(synonyms, str):
                synonyms = [synonyms]

            # string score
            if label_lower == term_lower or any(str(s).strip().lower() == term_lower for s in synonyms):
                string_score = 1.0
            else:
                string_score = SequenceMatcher(None, term_lower, label_lower).ratio()

            # keyword score
            if src_keywords:
                target_text = " ".join([label] + [str(s) for s in synonyms])
                tgt_keywords = self._extract_keywords(target_text)
                keyword_score = (len(src_keywords & tgt_keywords) / len(src_keywords)) if tgt_keywords else 0.0
            else:
                keyword_score = 0.0

            alpha = 0.7
            combined = alpha * string_score + (1.0 - alpha) * keyword_score
            iri = d.get("iri") or ""
            curie = iri_to_curie(iri)

            if combined > best[3]:
                best = (curie, iri, label, float(combined))

        if best[3] < float(min_combined):
            return ("", "", "", 0.0)

        if self.OLS_SLEEP_SECONDS and self.OLS_SLEEP_SECONDS > 0:
            time.sleep(self.OLS_SLEEP_SECONDS)

        return best

    def _bulk_ols_best(
        self,
        terms: List[str],
        ontology: str,
        ols_max_workers: int = 4,
        min_combined: float = 0.35,
    ) -> Dict[str, Tuple[str, str, str, float]]:
        """
        Map unique term strings -> (curie, iri, label, score).
        """
        terms = [t for t in terms if t]
        if not terms:
            return {}

        uniq_terms = sorted(set(terms))
        out: Dict[str, Tuple[str, str, str, float]] = {}

        if ols_max_workers <= 1:
            for t in uniq_terms:
                out[t] = self._ols_best_match(t, ontology=ontology, min_combined=min_combined)
            return out

        from concurrent.futures import ThreadPoolExecutor, as_completed
        with ThreadPoolExecutor(max_workers=ols_max_workers) as ex:
            futs = {
                ex.submit(self._ols_best_match, t, ontology, min_combined): t
                for t in uniq_terms
            }
            for fut in as_completed(futs):
                t = futs[fut]
                try:
                    out[t] = fut.result()
                except Exception:  # pragma: no cover
                    out[t] = ("", "", "", 0.0)

        return out

    def enrich_experimental_contexts(
        self,
        input_name: str = "ExperimentalContext_Properties_Processed.csv",
        output_name: str = "ExperimentalContext_Properties_WithOntologies.csv",
        ols_max_workers: int = 4,
    ) -> Optional[pd.DataFrame]:
        """
        Optional enrichment for ExperimentalContext nodes (FigureB):
          - Organism/species -> NCBITaxon
          - Cell line -> Cellosaurus
          - Tissue -> UBERON
          - Unit (if present) -> UO
        """
        input_path = f"{self.data_dir}/{input_name}"
        try:
            df = pd.read_csv(input_path, low_memory=False)
        except FileNotFoundError:
            logger.info("ExperimentalContext file not found: %s (skipping)", input_path)
            return None

        org_col = self._detect_first_column(df, ["Organism", "Species", "Taxon", "TaxonName", "OrganismName"])
        cell_col = self._detect_first_column(df, ["CellLine", "Cell_Line", "CellLineName"])
        tissue_col = self._detect_first_column(df, ["Tissue", "TissueName", "AnatomicalSite"])
        unit_col = self._detect_first_column(df, ["Unit", "UnitLabel", "UnitName"])

        if org_col:
            org_terms = df[org_col].fillna("").astype(str).str.strip().tolist()
            org_map = self._bulk_ols_best(org_terms, "ncbitaxon", ols_max_workers=ols_max_workers, min_combined=0.35)
            df["NCBITaxon_ID"] = [org_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in org_terms]
            df["NCBITaxon_IRI"] = [org_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in org_terms]
            df["NCBITaxon_Label"] = [org_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in org_terms]
            df["NCBITaxon_Score"] = [org_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in org_terms]
        else:
            df["NCBITaxon_ID"] = ""
            df["NCBITaxon_IRI"] = ""
            df["NCBITaxon_Label"] = ""
            df["NCBITaxon_Score"] = 0.0

        if cell_col:
            cell_terms = df[cell_col].fillna("").astype(str).str.strip().tolist()
            cell_map = self._bulk_ols_best(cell_terms, "cellosaurus", ols_max_workers=ols_max_workers, min_combined=0.35)
            df["Cellosaurus_ID"] = [cell_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in cell_terms]
            df["Cellosaurus_IRI"] = [cell_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in cell_terms]
            df["Cellosaurus_Label"] = [cell_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in cell_terms]
            df["Cellosaurus_Score"] = [cell_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in cell_terms]
        else:
            df["Cellosaurus_ID"] = ""
            df["Cellosaurus_IRI"] = ""
            df["Cellosaurus_Label"] = ""
            df["Cellosaurus_Score"] = 0.0

        if tissue_col:
            tissue_terms = df[tissue_col].fillna("").astype(str).str.strip().tolist()
            tissue_map = self._bulk_ols_best(tissue_terms, "uberon", ols_max_workers=ols_max_workers, min_combined=0.35)
            df["UBERON_ID"] = [tissue_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in tissue_terms]
            df["UBERON_IRI"] = [tissue_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in tissue_terms]
            df["UBERON_Label"] = [tissue_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in tissue_terms]
            df["UBERON_Score"] = [tissue_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in tissue_terms]
        else:
            df["UBERON_ID"] = ""
            df["UBERON_IRI"] = ""
            df["UBERON_Label"] = ""
            df["UBERON_Score"] = 0.0

        if unit_col:
            unit_terms = df[unit_col].fillna("").astype(str).str.strip().tolist()
            uo_map = self._bulk_ols_best(unit_terms, "uo", ols_max_workers=ols_max_workers, min_combined=0.30)
            df["UO_ID"] = [uo_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in unit_terms]
            df["UO_IRI"] = [uo_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in unit_terms]
            df["UO_Label"] = [uo_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in unit_terms]
            df["UO_Score"] = [uo_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in unit_terms]
            # Best-effort UCUM: if the value already looks like a UCUM code, keep it as-is
            df["UCUM_Code"] = [t if re.search(r"[A-Za-z].*", t) else "" for t in unit_terms]
        else:
            df["UO_ID"] = ""
            df["UO_IRI"] = ""
            df["UO_Label"] = ""
            df["UO_Score"] = 0.0
            df["UCUM_Code"] = ""

        df["OntologyEnrichedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.ExperimentalContextOLS/1.0"
        df["OntologySources"] = "OLS(ncbitaxon;cellosaurus;uberon;uo)"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved ExperimentalContext enriched table to %s", out_path)
        
        # Write a run summary JSON (mirrors the 4-core-node summaries)
        try:
            key_fields = ["NCBITaxon_ID", "Cellosaurus_ID", "UBERON_ID", "UO_ID"]
            def _nonempty(col: str) -> pd.Series:
                if col not in df.columns:
                    return pd.Series([False] * len(df))
                return df[col].fillna("").astype(str).str.strip().ne("")
            hit_any_mask = _nonempty("NCBITaxon_ID") | _nonempty("Cellosaurus_ID") | _nonempty("UBERON_ID") | _nonempty("UO_ID")
            summary = {
                "total_nodes_processed": int(len(df)),
                "total_successful_hits_any": int(hit_any_mask.sum()),
                "total_missing_entries": int((~hit_any_mask).sum()),
                "ontology_counts": {k: int(_nonempty(k).sum()) for k in key_fields},
                "ontology_sources": "OLS(ncbitaxon;cellosaurus;uberon;uo)",
                "generated_at_utc": datetime.datetime.utcnow().isoformat() + "Z",
                "input_file": input_name,
                "output_file": output_name,
            }
            summary_path = f"{self.data_dir}/ExperimentalContext_Ontology_RunSummary.json"
            with open(summary_path, "w", encoding="utf-8") as f:
                json.dump(summary, f, indent=2)
            logger.info("Saved ExperimentalContext ontology run summary to %s", summary_path)
        except Exception as exc:  # pragma: no cover
            logger.warning("Could not write ExperimentalContext run summary: %s", exc)

        return df

    def enrich_assay_endpoints(
        self,
        input_name: str = "AssayEndpoint_Properties_Processed.csv",
        output_name: str = "AssayEndpoint_Properties_WithOntologies.csv",
        endpoint_col: Optional[str] = None,
        ols_max_workers: int = 4,
    ) -> Optional[pd.DataFrame]:
        """
        Optional enrichment for AssayEndpoint nodes (FigureB):
          - endpoint_name -> BAO/EFO
          - unit -> UO (+ best-effort UCUM)
        """
        input_path = f"{self.data_dir}/{input_name}"
        try:
            df = pd.read_csv(input_path, low_memory=False)
        except FileNotFoundError:
            logger.info("AssayEndpoint file not found: %s (skipping)", input_path)
            return None

        if endpoint_col is None:
            endpoint_col = self._detect_first_column(df, ["endpoint_name", "EndpointName", "Endpoint", "ReadoutType"])
        if endpoint_col is None:
            endpoint_col = df.columns[0]

        endpoint_terms = df[endpoint_col].fillna("").astype(str).str.strip().tolist()
        bao_map = self._bulk_ols_best(endpoint_terms, "bao", ols_max_workers=ols_max_workers, min_combined=0.35)
        efo_map = self._bulk_ols_best(endpoint_terms, "efo", ols_max_workers=ols_max_workers, min_combined=0.35)

        df["BAO_ID"] = [bao_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in endpoint_terms]
        df["BAO_IRI"] = [bao_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in endpoint_terms]
        df["BAO_Label"] = [bao_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in endpoint_terms]
        df["BAO_Score"] = [bao_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in endpoint_terms]

        df["EFO_ID"] = [efo_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in endpoint_terms]
        df["EFO_IRI"] = [efo_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in endpoint_terms]
        df["EFO_Label"] = [efo_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in endpoint_terms]
        df["EFO_Score"] = [efo_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in endpoint_terms]

        unit_col = self._detect_first_column(df, ["unit", "Unit", "UnitLabel", "UnitName"])
        if unit_col:
            unit_terms = df[unit_col].fillna("").astype(str).str.strip().tolist()
            uo_map = self._bulk_ols_best(unit_terms, "uo", ols_max_workers=ols_max_workers, min_combined=0.30)
            df["UO_ID"] = [uo_map.get(t, ("", "", "", 0.0))[0] if t else "" for t in unit_terms]
            df["UO_IRI"] = [uo_map.get(t, ("", "", "", 0.0))[1] if t else "" for t in unit_terms]
            df["UO_Label"] = [uo_map.get(t, ("", "", "", 0.0))[2] if t else "" for t in unit_terms]
            df["UO_Score"] = [uo_map.get(t, ("", "", "", 0.0))[3] if t else 0.0 for t in unit_terms]
            df["UCUM_Code"] = [t if re.search(r"[A-Za-z].*", t) else "" for t in unit_terms]
        else:
            df["UO_ID"] = ""
            df["UO_IRI"] = ""
            df["UO_Label"] = ""
            df["UO_Score"] = 0.0
            df["UCUM_Code"] = ""

        df["OntologyEnrichedAt"] = datetime.datetime.utcnow().isoformat() + "Z"
        df["OntologyEnricherVersion"] = "ChemGraphBuilder.AssayEndpointOLS/1.0"
        df["OntologySources"] = "OLS(bao;efo;uo)"

        out_path = f"{self.data_dir}/{output_name}"
        df.to_csv(out_path, index=False)
        logger.info("Saved AssayEndpoint enriched table to %s", out_path)
        
        # Write a run summary JSON (mirrors the 4-core-node summaries)
        try:
            key_fields = ["BAO_ID", "EFO_ID", "UO_ID"]
            def _nonempty(col: str) -> pd.Series:
                if col not in df.columns:
                    return pd.Series([False] * len(df))
                return df[col].fillna("").astype(str).str.strip().ne("")
            hit_any_mask = _nonempty("BAO_ID") | _nonempty("EFO_ID") | _nonempty("UO_ID")
            summary = {
                "total_nodes_processed": int(len(df)),
                "total_successful_hits_any": int(hit_any_mask.sum()),
                "total_missing_entries": int((~hit_any_mask).sum()),
                "ontology_counts": {k: int(_nonempty(k).sum()) for k in key_fields},
                "ontology_sources": "OLS(bao;efo;uo)",
                "generated_at_utc": datetime.datetime.utcnow().isoformat() + "Z",
                "input_file": input_name,
                "output_file": output_name,
            }
            summary_path = f"{self.data_dir}/AssayEndpoint_Ontology_RunSummary.json"
            with open(summary_path, "w", encoding="utf-8") as f:
                json.dump(summary, f, indent=2)
            logger.info("Saved AssayEndpoint ontology run summary to %s", summary_path)
        except Exception as exc:  # pragma: no cover
            logger.warning("Could not write AssayEndpoint run summary: %s", exc)

        return df
