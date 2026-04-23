#!/usr/bin/env python3
"""Extract yeast polymerase transcription factor domain annotations from SGD.

Pipeline:
1. Parse the Excel workbook without third-party dependencies.
2. Read the combined protein list (P:S) and map each protein back to Pol I/II/III.
3. Query SGD for raw domain annotations.
4. Query InterPro per protein to harmonize source-specific domain IDs into InterPro entries.
5. Write raw, harmonized, source-summary, and matrix outputs.

Default input workbook:
    ./data/input/polymerase_tf_sequences.xlsx

Alternative CSV input:
    ./data/input/proteins.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen
from zipfile import ZipFile
import xml.etree.ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_XLSX = REPO_ROOT / "data" / "input" / "polymerase_tf_sequences.xlsx"
DEFAULT_CSV = REPO_ROOT / "data" / "input" / "proteins.csv"
DEFAULT_OUTDIR = REPO_ROOT / "data" / "output"
SGD_BASE = "https://www.yeastgenome.org/backend/locus"
INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot"
NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

CATEGORY_BLOCKS = {
    "Pol I": ("A", "B", "C", "D"),
    "Pol II": ("F", "G", "H", "I"),
    "Pol III": ("K", "L", "M", "N"),
}
COMBINED_BLOCK = ("P", "Q", "R", "S")
INTERPRO_MEMBER_TO_SOURCE = {
    "cathgene3d": "Gene3D",
    "cdd": "CDD",
    "mobidblite": "MobiDBLite",
    "ncbifam": "TIGRFAM",
    "panther": "PANTHER",
    "pfam": "Pfam",
    "pirsf": "PIRSF",
    "prints": "PRINTS",
    "profile": "PROSITE",
    "prosite": "PROSITE",
    "smart": "SMART",
    "ssf": "SUPERFAMILY",
}

CSV_FIELD_ALIASES = {
    "category": ("category", "group", "polymerase", "pol", "class"),
    "protein": ("protein", "protein_name", "name", "gene", "standard_name"),
    "y_name": ("y_name", "y name", "yname", "systematic_name", "systematic name", "locus", "sgd_locus"),
    "accession": ("accession", "accession_id", "uniprot", "uniprot_id"),
    "fasta_header": ("fasta_header", "fasta header", "header"),
    "amino_acid_sequence": (
        "amino_acid_sequence",
        "amino acid sequence",
        "sequence",
        "protein_sequence",
        "aa_sequence",
    ),
    "raw_fasta": ("raw_fasta", "raw fasta", "fasta"),
}


def col_to_num(col: str) -> int:
    n = 0
    for ch in col:
        if ch.isalpha():
            n = n * 26 + ord(ch.upper()) - 64
    return n


def clean_text(value: str | None) -> str:
    if value is None:
        return ""
    return str(value).replace("\xa0", " ").strip()


def join_sorted(values: Iterable[str], sep: str = "; ") -> str:
    cleaned = sorted({clean_text(v) for v in values if clean_text(v)})
    return sep.join(cleaned)


def parse_fasta_cell(raw: str) -> Tuple[str, str]:
    raw = clean_text(raw)
    if not raw:
        return "", ""
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    header = lines[0] if lines and lines[0].startswith(">") else ""
    seq_lines = lines[1:] if header else lines
    sequence = "".join(seq_lines).replace(" ", "")
    return header, sequence


def normalize_csv_header(value: str | None) -> str:
    return "".join(ch for ch in clean_text(value).lower() if ch.isalnum())


def csv_value(row: dict, field: str) -> str:
    for alias in CSV_FIELD_ALIASES[field]:
        value = row.get(normalize_csv_header(alias), "")
        if clean_text(value):
            return clean_text(value)
    return ""


def read_csv_entries(path: Path) -> List[dict]:
    entries: List[dict] = []
    seen_y_names: set[str] = set()

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SystemExit(f"CSV input has no header row: {path}")

        normalized_fields = {normalize_csv_header(field): field for field in reader.fieldnames}
        has_lookup_field = any(
            normalize_csv_header(alias) in normalized_fields
            for field in ("y_name", "protein")
            for alias in CSV_FIELD_ALIASES[field]
        )
        if not has_lookup_field:
            expected = ", ".join((*CSV_FIELD_ALIASES["y_name"], *CSV_FIELD_ALIASES["protein"]))
            raise SystemExit(f"CSV input must include an SGD lookup column such as: {expected}")

        for row_num, raw_row in enumerate(reader, start=2):
            row = {normalize_csv_header(key): value for key, value in raw_row.items() if key is not None}
            y_name = csv_value(row, "y_name") or csv_value(row, "protein")
            if not y_name:
                continue
            if y_name in seen_y_names:
                raise SystemExit(f"Duplicate y_name in CSV input at row {row_num}: {y_name}")
            seen_y_names.add(y_name)

            raw_fasta = csv_value(row, "raw_fasta")
            fasta_header, fasta_sequence = parse_fasta_cell(raw_fasta)
            header = csv_value(row, "fasta_header") or fasta_header
            sequence = csv_value(row, "amino_acid_sequence") or fasta_sequence
            protein = csv_value(row, "protein") or y_name

            entries.append(
                {
                    "category": csv_value(row, "category") or "Unknown",
                    "protein": protein,
                    "y_name": y_name,
                    "accession": csv_value(row, "accession"),
                    "fasta_header": header,
                    "amino_acid_sequence": sequence.replace(" ", "").replace("\n", ""),
                    "raw_fasta": raw_fasta,
                }
            )

    if not entries:
        raise SystemExit(f"No protein rows found in CSV input: {path}")
    return entries


class XlsxSheet:
    def __init__(self, path: Path):
        self.path = path
        self.rows = self._load_rows()

    def _load_rows(self) -> Dict[int, Dict[int, str]]:
        rows: Dict[int, Dict[int, str]] = defaultdict(dict)
        with ZipFile(self.path) as zf:
            shared_strings: List[str] = []
            if "xl/sharedStrings.xml" in zf.namelist():
                root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
                for si in root.findall("a:si", NS):
                    parts = [t.text or "" for t in si.iterfind(".//a:t", NS)]
                    shared_strings.append("".join(parts))

            sheet = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
            for cell in sheet.iterfind(".//a:sheetData/a:row/a:c", NS):
                ref = cell.attrib.get("r", "")
                letters = "".join(ch for ch in ref if ch.isalpha())
                digits = "".join(ch for ch in ref if ch.isdigit())
                if not letters or not digits:
                    continue
                row_num = int(digits)
                col_num = col_to_num(letters)
                value_node = cell.find("a:v", NS)
                if value_node is None:
                    continue
                value = value_node.text or ""
                if cell.attrib.get("t") == "s":
                    value = shared_strings[int(value)]
                rows[row_num][col_num] = value
        return rows

    def get(self, col: str, row: int) -> str:
        return clean_text(self.rows.get(row, {}).get(col_to_num(col)))

    @property
    def max_row(self) -> int:
        return max(self.rows) if self.rows else 0


def extract_category_entries(sheet: XlsxSheet) -> Tuple[List[dict], Dict[str, str]]:
    proteins: List[dict] = []
    category_by_yname: Dict[str, str] = {}
    for category, (c_protein, c_yname, c_acc, c_fasta) in CATEGORY_BLOCKS.items():
        for row in range(3, sheet.max_row + 1):
            protein = sheet.get(c_protein, row)
            y_name = sheet.get(c_yname, row)
            accession = sheet.get(c_acc, row)
            raw_fasta = sheet.get(c_fasta, row)
            if not protein:
                continue
            header, sequence = parse_fasta_cell(raw_fasta)
            entry = {
                "category": category,
                "protein": protein,
                "y_name": y_name,
                "accession": accession,
                "fasta_header": header,
                "amino_acid_sequence": sequence,
                "raw_fasta": raw_fasta,
            }
            proteins.append(entry)
            if y_name:
                category_by_yname[y_name] = category
    return proteins, category_by_yname


def extract_combined_entries(sheet: XlsxSheet, category_by_yname: Dict[str, str]) -> List[dict]:
    combined: List[dict] = []
    c_protein, c_yname, c_acc, c_fasta = COMBINED_BLOCK
    for row in range(3, sheet.max_row + 1):
        protein = sheet.get(c_protein, row)
        y_name = sheet.get(c_yname, row)
        accession = sheet.get(c_acc, row)
        raw_fasta = sheet.get(c_fasta, row)
        if not protein:
            continue
        header, sequence = parse_fasta_cell(raw_fasta)
        combined.append(
            {
                "category": category_by_yname.get(y_name, "Unknown"),
                "protein": protein,
                "y_name": y_name,
                "accession": accession,
                "fasta_header": header,
                "amino_acid_sequence": sequence,
                "raw_fasta": raw_fasta,
            }
        )
    return combined


def fetch_json(url: str, timeout: float = 30.0) -> object | None:
    req = Request(url, headers={"User-Agent": "yeast-motif-pipeline/1.0"})
    with urlopen(req, timeout=timeout) as response:
        body = response.read().decode("utf-8")
        if not body.strip():
            return None
        return json.loads(body)


def cached_json(url: str, cache_path: Path, force_refresh: bool, timeout: float) -> object | None:
    if cache_path.exists() and not force_refresh:
        return json.loads(cache_path.read_text())
    data = fetch_json(url, timeout=timeout)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(data, indent=2, sort_keys=True))
    return data


def fetch_sgd_records(
    proteins: List[dict],
    outdir: Path,
    sleep_seconds: float,
    force_refresh: bool,
    timeout: float,
) -> Tuple[List[dict], List[dict], List[dict]]:
    cache_dir = outdir / "cache" / "sgd"
    protein_rows: List[dict] = []
    domain_rows: List[dict] = []
    failures: List[dict] = []

    for idx, protein in enumerate(proteins, start=1):
        y_name = protein["y_name"]
        if not y_name:
            failures.append({"y_name": "", "protein": protein["protein"], "error": "Missing Y Name"})
            continue

        encoded = quote(y_name, safe="")
        locus_url = f"{SGD_BASE}/{encoded}"
        domains_url = f"{SGD_BASE}/{encoded}/protein_domain_details"
        locus_cache = cache_dir / f"{y_name}.locus.json"
        domains_cache = cache_dir / f"{y_name}.domains.json"

        try:
            locus = cached_json(locus_url, locus_cache, force_refresh, timeout)
            domains = cached_json(domains_url, domains_cache, force_refresh, timeout)
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
            failures.append({"y_name": y_name, "protein": protein["protein"], "error": str(exc)})
            continue

        protein_row = dict(protein)
        if isinstance(locus, dict):
            protein_row.update(
                {
                    "sgdid": clean_text(locus.get("sgdid")),
                    "standard_name": clean_text(locus.get("gene_name") or locus.get("display_name")),
                    "systematic_name": clean_text(locus.get("format_name")),
                    "uniprot_id": clean_text(locus.get("uniprot_id")),
                    "sgd_link": f"https://www.yeastgenome.org{clean_text(locus.get('link'))}" if clean_text(locus.get('link')) else "",
                }
            )
        protein_rows.append(protein_row)

        if isinstance(domains, list):
            for hit in domains:
                domain = hit.get("domain", {}) if isinstance(hit, dict) else {}
                source = hit.get("source", {}) if isinstance(hit, dict) else {}
                locus_part = hit.get("locus", {}) if isinstance(hit, dict) else {}
                domain_id = clean_text(domain.get("display_name"))
                source_name = clean_text(source.get("display_name") or source.get("format_name"))
                raw_domain_key = f"{source_name}::{domain_id}" if source_name and domain_id else ""
                domain_rows.append(
                    {
                        "category": protein["category"],
                        "protein": protein["protein"],
                        "y_name": protein["y_name"],
                        "accession": protein["accession"],
                        "sgdid": protein_row.get("sgdid", ""),
                        "standard_name": protein_row.get("standard_name", ""),
                        "systematic_name": protein_row.get("systematic_name", protein["y_name"]),
                        "uniprot_id": protein_row.get("uniprot_id", ""),
                        "domain_hit_id": hit.get("id", ""),
                        "source": source_name,
                        "source_format_name": clean_text(source.get("format_name")),
                        "raw_domain_id": domain_id,
                        "raw_domain_key": raw_domain_key,
                        "domain_description": clean_text(domain.get("description")),
                        "domain_link": f"https://www.yeastgenome.org{clean_text(domain.get('link'))}" if clean_text(domain.get('link')) else "",
                        "domain_gene_count": domain.get("count", ""),
                        "start": hit.get("start", ""),
                        "end": hit.get("end", ""),
                        "locus_display_name": clean_text(locus_part.get("display_name")),
                        "locus_format_name": clean_text(locus_part.get("format_name")),
                    }
                )

        if sleep_seconds and idx < len(proteins):
            time.sleep(sleep_seconds)

    return protein_rows, domain_rows, failures


def fetch_interpro_harmonization(
    protein_rows: List[dict],
    outdir: Path,
    sleep_seconds: float,
    force_refresh: bool,
    timeout: float,
) -> Tuple[Dict[str, Dict[str, List[dict]]], List[dict]]:
    cache_dir = outdir / "cache" / "interpro"
    mappings: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: defaultdict(list))
    failures: List[dict] = []

    for idx, protein in enumerate(protein_rows, start=1):
        y_name = protein.get("y_name", "")
        uniprot_id = protein.get("uniprot_id", "")
        if not y_name or not uniprot_id:
            failures.append({"y_name": y_name, "protein": protein.get("protein", ""), "error": "Missing UniProt ID"})
            continue

        url = f"{INTERPRO_BASE}/{quote(uniprot_id, safe='')}"
        cache_path = cache_dir / f"{y_name}.interpro.json"
        try:
            data = cached_json(url, cache_path, force_refresh, timeout)
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
            failures.append({"y_name": y_name, "protein": protein.get("protein", ""), "error": str(exc)})
            continue

        if not data:
            if sleep_seconds and idx < len(protein_rows):
                time.sleep(sleep_seconds)
            continue

        results = data.get("results", []) if isinstance(data, dict) else []
        for result in results:
            metadata = result.get("metadata", {})
            interpro_accession = clean_text(metadata.get("accession"))
            interpro_name = clean_text(metadata.get("name"))
            interpro_type = clean_text(metadata.get("type"))
            member_databases = metadata.get("member_databases") or {}
            for member_db, accession_map in member_databases.items():
                sgd_source = INTERPRO_MEMBER_TO_SOURCE.get(member_db.lower())
                if not sgd_source or not isinstance(accession_map, dict):
                    continue
                for member_accession, member_name in accession_map.items():
                    raw_key = f"{sgd_source}::{clean_text(member_accession)}"
                    mappings[y_name][raw_key].append(
                        {
                            "interpro_accession": interpro_accession,
                            "interpro_name": interpro_name,
                            "interpro_type": interpro_type,
                            "member_database": member_db,
                            "member_accession": clean_text(member_accession),
                            "member_name": clean_text(member_name),
                        }
                    )

        if sleep_seconds and idx < len(protein_rows):
            time.sleep(sleep_seconds)

    return mappings, failures


def annotate_harmonization(domain_rows: List[dict], interpro_mappings: Dict[str, Dict[str, List[dict]]]) -> List[dict]:
    annotated: List[dict] = []
    for row in domain_rows:
        raw_key = row.get("raw_domain_key", "")
        y_name = row.get("y_name", "")
        candidates = interpro_mappings.get(y_name, {}).get(raw_key, [])
        interpro_accessions = sorted({c["interpro_accession"] for c in candidates if c.get("interpro_accession")})

        if len(interpro_accessions) == 1:
            harmonization_status = "mapped"
            harmonized_accession = interpro_accessions[0]
            chosen = next(c for c in candidates if c.get("interpro_accession") == harmonized_accession)
            harmonized_name = chosen.get("interpro_name", "")
            harmonized_type = chosen.get("interpro_type", "")
            harmonized_key = f"InterPro::{harmonized_accession}"
        elif len(interpro_accessions) > 1:
            harmonization_status = "ambiguous"
            harmonized_accession = join_sorted(interpro_accessions, sep=";")
            harmonized_name = join_sorted(c.get("interpro_name", "") for c in candidates)
            harmonized_type = join_sorted(c.get("interpro_type", "") for c in candidates)
            harmonized_key = row.get("raw_domain_key", "")
        else:
            harmonization_status = "unmapped"
            harmonized_accession = ""
            harmonized_name = ""
            harmonized_type = ""
            harmonized_key = row.get("raw_domain_key", "")

        annotated_row = dict(row)
        annotated_row.update(
            {
                "interpro_candidate_count": len(interpro_accessions),
                "interpro_candidates": join_sorted(interpro_accessions, sep=";"),
                "harmonization_status": harmonization_status,
                "harmonized_accession": harmonized_accession,
                "harmonized_name": harmonized_name,
                "harmonized_type": harmonized_type,
                "harmonized_key": harmonized_key,
                "interpro_member_sources": join_sorted(c.get("member_database", "") for c in candidates),
                "interpro_member_accessions": join_sorted(c.get("member_accession", "") for c in candidates),
            }
        )
        annotated.append(annotated_row)
    return annotated


def collapse_harmonized_hits(annotated_rows: List[dict]) -> List[dict]:
    grouped: Dict[Tuple[str, str], dict] = {}
    for row in annotated_rows:
        key = (row.get("y_name", ""), row.get("harmonized_key", ""))
        group = grouped.setdefault(
            key,
            {
                "category": row.get("category", ""),
                "protein": row.get("protein", ""),
                "y_name": row.get("y_name", ""),
                "accession": row.get("accession", ""),
                "sgdid": row.get("sgdid", ""),
                "standard_name": row.get("standard_name", ""),
                "systematic_name": row.get("systematic_name", ""),
                "uniprot_id": row.get("uniprot_id", ""),
                "harmonization_status": row.get("harmonization_status", ""),
                "harmonized_key": row.get("harmonized_key", ""),
                "harmonized_accession": row.get("harmonized_accession", ""),
                "harmonized_name": row.get("harmonized_name", ""),
                "harmonized_type": row.get("harmonized_type", ""),
                "support_hit_count": 0,
                "support_sources": set(),
                "support_raw_domain_ids": set(),
                "support_raw_domain_keys": set(),
                "support_domain_descriptions": set(),
                "support_ranges": set(),
            },
        )
        group["support_hit_count"] += 1
        group["support_sources"].add(row.get("source", ""))
        group["support_raw_domain_ids"].add(row.get("raw_domain_id", ""))
        group["support_raw_domain_keys"].add(row.get("raw_domain_key", ""))
        group["support_domain_descriptions"].add(row.get("domain_description", ""))
        start = clean_text(row.get("start", ""))
        end = clean_text(row.get("end", ""))
        if start or end:
            group["support_ranges"].add(f"{start}-{end}")

    collapsed: List[dict] = []
    for group in grouped.values():
        collapsed.append(
            {
                **{k: v for k, v in group.items() if not isinstance(v, set)},
                "support_sources": join_sorted(group["support_sources"]),
                "support_raw_domain_ids": join_sorted(group["support_raw_domain_ids"]),
                "support_raw_domain_keys": join_sorted(group["support_raw_domain_keys"]),
                "support_domain_descriptions": join_sorted(group["support_domain_descriptions"]),
                "support_ranges": join_sorted(group["support_ranges"]),
            }
        )
    collapsed.sort(key=lambda row: (row.get("category", ""), row.get("protein", ""), row.get("harmonized_key", "")))
    return collapsed


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_domain_matrices(
    domain_rows: List[dict], proteins: List[dict], key_field: str, min_count: int, count_field: str | None = None
) -> Tuple[List[dict], List[dict], List[dict], List[dict], List[str], List[str]]:
    protein_order = [protein["y_name"] for protein in proteins]
    labels = {protein["y_name"]: protein for protein in proteins}

    domain_counts_by_protein: Dict[str, Counter] = defaultdict(Counter)
    domain_presence_counts: Counter = Counter()
    for row in domain_rows:
        domain_key = clean_text(row.get(key_field, ""))
        y_name = clean_text(row.get("y_name", ""))
        if not domain_key or not y_name:
            continue
        increment = row.get(count_field, 1) if count_field else 1
        try:
            increment = int(increment)
        except (TypeError, ValueError):
            increment = 1
        domain_counts_by_protein[y_name][domain_key] += increment

    for counter in domain_counts_by_protein.values():
        for domain_key in counter:
            domain_presence_counts[domain_key] += 1

    all_domains = sorted(domain_presence_counts)
    filtered_domains = [domain for domain in all_domains if domain_presence_counts[domain] >= min_count]

    def matrix_rows(domain_list: List[str], mode: str) -> List[dict]:
        rows: List[dict] = []
        for y_name in protein_order:
            protein = labels[y_name]
            row = {
                "category": protein["category"],
                "protein": protein["protein"],
                "y_name": protein["y_name"],
            }
            counts = domain_counts_by_protein.get(y_name, Counter())
            for domain in domain_list:
                value = counts.get(domain, 0)
                row[domain] = 1 if mode == "binary" and value > 0 else value
            rows.append(row)
        return rows

    return (
        matrix_rows(all_domains, mode="binary"),
        matrix_rows(filtered_domains, mode="binary"),
        matrix_rows(all_domains, mode="count"),
        matrix_rows(filtered_domains, mode="count"),
        all_domains,
        filtered_domains,
    )


def build_source_summary(annotated_rows: List[dict]) -> List[dict]:
    grouped: Dict[str, dict] = {}
    for row in annotated_rows:
        source = row.get("source", "") or "Unknown"
        entry = grouped.setdefault(
            source,
            {
                "source": source,
                "domain_hits": 0,
                "unique_raw_domain_ids": set(),
                "unique_raw_domain_keys": set(),
                "mapped_hits": 0,
                "mapped_raw_domain_keys": set(),
                "mapped_harmonized_keys": set(),
                "unmapped_hits": 0,
                "unmapped_raw_domain_keys": set(),
                "ambiguous_hits": 0,
                "ambiguous_raw_domain_keys": set(),
            },
        )
        entry["domain_hits"] += 1
        entry["unique_raw_domain_ids"].add(row.get("raw_domain_id", ""))
        entry["unique_raw_domain_keys"].add(row.get("raw_domain_key", ""))

        status = row.get("harmonization_status", "")
        if status == "mapped":
            entry["mapped_hits"] += 1
            entry["mapped_raw_domain_keys"].add(row.get("raw_domain_key", ""))
            entry["mapped_harmonized_keys"].add(row.get("harmonized_key", ""))
        elif status == "ambiguous":
            entry["ambiguous_hits"] += 1
            entry["ambiguous_raw_domain_keys"].add(row.get("raw_domain_key", ""))
        else:
            entry["unmapped_hits"] += 1
            entry["unmapped_raw_domain_keys"].add(row.get("raw_domain_key", ""))

    summary_rows: List[dict] = []
    for source, entry in sorted(grouped.items()):
        summary_rows.append(
            {
                "source": source,
                "domain_hits": entry["domain_hits"],
                "unique_raw_domain_ids": len({v for v in entry["unique_raw_domain_ids"] if v}),
                "unique_raw_domain_keys": len({v for v in entry["unique_raw_domain_keys"] if v}),
                "mapped_hits": entry["mapped_hits"],
                "mapped_unique_raw_domain_keys": len({v for v in entry["mapped_raw_domain_keys"] if v}),
                "mapped_unique_harmonized_keys": len({v for v in entry["mapped_harmonized_keys"] if v}),
                "unmapped_hits": entry["unmapped_hits"],
                "unmapped_unique_raw_domain_keys": len({v for v in entry["unmapped_raw_domain_keys"] if v}),
                "ambiguous_hits": entry["ambiguous_hits"],
                "ambiguous_unique_raw_domain_keys": len({v for v in entry["ambiguous_raw_domain_keys"] if v}),
            }
        )
    return summary_rows


def summarize_run(protein_count: int, domain_rows: List[dict], key_field: str, min_count: int, failures: int) -> dict:
    domains_by_protein: Dict[str, set] = defaultdict(set)
    for row in domain_rows:
        key = clean_text(row.get(key_field, ""))
        y_name = clean_text(row.get("y_name", ""))
        if key and y_name:
            domains_by_protein[y_name].add(key)
    domain_counts = Counter()
    for keys in domains_by_protein.values():
        for key in keys:
            domain_counts[key] += 1
    return {
        "protein_count": protein_count,
        "domain_hit_count": len(domain_rows),
        "unique_domains": len(domain_counts),
        "unique_domains_min_count": sum(1 for count in domain_counts.values() if count >= min_count),
        "min_proteins_per_domain": min_count,
        "failure_count": failures,
    }


def validate_workbook_counts(category_entries: List[dict], combined_entries: List[dict]) -> None:
    counts = Counter(entry["category"] for entry in category_entries)
    expected = {"Pol I": 13, "Pol II": 26, "Pol III": 9}
    if dict(counts) != expected:
        raise SystemExit(f"Unexpected category counts: {dict(counts)} != {expected}")
    if len(combined_entries) != 48:
        raise SystemExit(f"Expected 48 combined proteins, found {len(combined_entries)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="Path to the input .xlsx workbook")
    input_group.add_argument(
        "--csv",
        type=Path,
        help=(
            "Path to a CSV protein table. Preferred lookup column: y_name "
            "(or Y Name/systematic_name/locus); protein/name/standard_name also work as lookup columns. "
            "Optional columns: category, protein, accession, fasta_header, amino_acid_sequence, raw_fasta."
        ),
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Output directory for CSV and cache files")
    parser.add_argument("--sleep-seconds", type=float, default=0.15, help="Delay between API requests")
    parser.add_argument("--timeout", type=float, default=30.0, help="HTTP timeout in seconds")
    parser.add_argument("--min-proteins-per-domain", type=int, default=2, help="Minimum protein count for filtered matrices")
    parser.add_argument("--force-refresh", action="store_true", help="Ignore cached API responses and fetch again")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.csv or args.xlsx
    if not input_path.exists():
        input_type = "CSV input" if args.csv else "Input workbook"
        print(f"{input_type} not found: {input_path}", file=sys.stderr)
        return 1

    if args.csv:
        combined_entries = read_csv_entries(args.csv)
    else:
        sheet = XlsxSheet(args.xlsx)
        category_entries, category_by_yname = extract_category_entries(sheet)
        combined_entries = extract_combined_entries(sheet, category_by_yname)
        validate_workbook_counts(category_entries, combined_entries)

    protein_rows, raw_domain_rows, sgd_failures = fetch_sgd_records(
        proteins=combined_entries,
        outdir=args.outdir,
        sleep_seconds=args.sleep_seconds,
        force_refresh=args.force_refresh,
        timeout=args.timeout,
    )
    interpro_mappings, interpro_failures = fetch_interpro_harmonization(
        protein_rows=protein_rows,
        outdir=args.outdir,
        sleep_seconds=args.sleep_seconds,
        force_refresh=args.force_refresh,
        timeout=args.timeout,
    )
    annotated_rows = annotate_harmonization(raw_domain_rows, interpro_mappings)
    collapsed_harmonized_rows = collapse_harmonized_hits(annotated_rows)
    source_summary_rows = build_source_summary(annotated_rows)

    (
        raw_all_matrix_rows,
        raw_filtered_matrix_rows,
        raw_all_count_matrix_rows,
        raw_filtered_count_matrix_rows,
        raw_all_domains,
        raw_filtered_domains,
    ) = build_domain_matrices(
        annotated_rows, combined_entries, key_field="raw_domain_key", min_count=args.min_proteins_per_domain
    )
    (
        harm_all_matrix_rows,
        harm_filtered_matrix_rows,
        harm_all_count_matrix_rows,
        harm_filtered_count_matrix_rows,
        harm_all_domains,
        harm_filtered_domains,
    ) = build_domain_matrices(
        collapsed_harmonized_rows,
        combined_entries,
        key_field="harmonized_key",
        min_count=args.min_proteins_per_domain,
        count_field="support_hit_count",
    )

    protein_fields = [
        "category",
        "protein",
        "y_name",
        "accession",
        "sgdid",
        "standard_name",
        "systematic_name",
        "uniprot_id",
        "sgd_link",
        "fasta_header",
        "amino_acid_sequence",
        "raw_fasta",
    ]
    raw_domain_fields = [
        "category",
        "protein",
        "y_name",
        "accession",
        "sgdid",
        "standard_name",
        "systematic_name",
        "uniprot_id",
        "domain_hit_id",
        "source",
        "source_format_name",
        "raw_domain_id",
        "raw_domain_key",
        "domain_description",
        "domain_link",
        "domain_gene_count",
        "start",
        "end",
        "locus_display_name",
        "locus_format_name",
        "interpro_candidate_count",
        "interpro_candidates",
        "harmonization_status",
        "harmonized_accession",
        "harmonized_name",
        "harmonized_type",
        "harmonized_key",
        "interpro_member_sources",
        "interpro_member_accessions",
    ]
    collapsed_fields = [
        "category",
        "protein",
        "y_name",
        "accession",
        "sgdid",
        "standard_name",
        "systematic_name",
        "uniprot_id",
        "harmonization_status",
        "harmonized_key",
        "harmonized_accession",
        "harmonized_name",
        "harmonized_type",
        "support_hit_count",
        "support_sources",
        "support_raw_domain_ids",
        "support_raw_domain_keys",
        "support_domain_descriptions",
        "support_ranges",
    ]
    source_summary_fields = [
        "source",
        "domain_hits",
        "unique_raw_domain_ids",
        "unique_raw_domain_keys",
        "mapped_hits",
        "mapped_unique_raw_domain_keys",
        "mapped_unique_harmonized_keys",
        "unmapped_hits",
        "unmapped_unique_raw_domain_keys",
        "ambiguous_hits",
        "ambiguous_unique_raw_domain_keys",
    ]
    failure_fields = ["y_name", "protein", "error"]

    outdir = args.outdir
    write_csv(outdir / "polymerase_tf_proteins.csv", protein_rows, protein_fields)
    write_csv(outdir / "polymerase_tf_domain_hits.csv", annotated_rows, raw_domain_fields)
    write_csv(outdir / "polymerase_tf_domain_hits_raw.csv", annotated_rows, raw_domain_fields)
    write_csv(outdir / "polymerase_tf_domain_hits_harmonized.csv", collapsed_harmonized_rows, collapsed_fields)
    write_csv(outdir / "polymerase_tf_domain_source_summary.csv", source_summary_rows, source_summary_fields)
    write_csv(outdir / "sgd_failures.csv", sgd_failures, failure_fields)
    write_csv(outdir / "interpro_failures.csv", interpro_failures, failure_fields)

    write_csv(
        outdir / "polymerase_tf_domain_matrix_all.csv",
        raw_all_matrix_rows,
        ["category", "protein", "y_name", *raw_all_domains],
    )
    write_csv(
        outdir / f"polymerase_tf_domain_matrix_min{args.min_proteins_per_domain}.csv",
        raw_filtered_matrix_rows,
        ["category", "protein", "y_name", *raw_filtered_domains],
    )
    write_csv(
        outdir / "polymerase_tf_domain_matrix_counts_all.csv",
        raw_all_count_matrix_rows,
        ["category", "protein", "y_name", *raw_all_domains],
    )
    write_csv(
        outdir / f"polymerase_tf_domain_matrix_counts_min{args.min_proteins_per_domain}.csv",
        raw_filtered_count_matrix_rows,
        ["category", "protein", "y_name", *raw_filtered_domains],
    )
    write_csv(
        outdir / "polymerase_tf_domain_matrix_harmonized_all.csv",
        harm_all_matrix_rows,
        ["category", "protein", "y_name", *harm_all_domains],
    )
    write_csv(
        outdir / f"polymerase_tf_domain_matrix_harmonized_min{args.min_proteins_per_domain}.csv",
        harm_filtered_matrix_rows,
        ["category", "protein", "y_name", *harm_filtered_domains],
    )
    write_csv(
        outdir / "polymerase_tf_domain_matrix_harmonized_counts_all.csv",
        harm_all_count_matrix_rows,
        ["category", "protein", "y_name", *harm_all_domains],
    )
    write_csv(
        outdir / f"polymerase_tf_domain_matrix_harmonized_counts_min{args.min_proteins_per_domain}.csv",
        harm_filtered_count_matrix_rows,
        ["category", "protein", "y_name", *harm_filtered_domains],
    )

    raw_summary = summarize_run(
        protein_count=len(combined_entries),
        domain_rows=annotated_rows,
        key_field="raw_domain_key",
        min_count=args.min_proteins_per_domain,
        failures=len(sgd_failures),
    )
    harmonized_summary = summarize_run(
        protein_count=len(combined_entries),
        domain_rows=collapsed_harmonized_rows,
        key_field="harmonized_key",
        min_count=args.min_proteins_per_domain,
        failures=len(sgd_failures) + len(interpro_failures),
    )
    summary = {
        "input_file": str(input_path),
        "input_format": "csv" if args.csv else "xlsx",
        "input_workbook": str(args.xlsx) if not args.csv else "",
        "input_csv": str(args.csv) if args.csv else "",
        "output_directory": str(outdir),
        "category_counts": Counter(entry["category"] for entry in combined_entries),
        "raw": raw_summary,
        "harmonized": harmonized_summary,
        "source_summary": source_summary_rows,
        "sgd_failure_count": len(sgd_failures),
        "interpro_failure_count": len(interpro_failures),
    }
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, default=dict))

    print("Done.")
    print("\nRaw run results")
    print(f"Proteins: {raw_summary['protein_count']}")
    print(f"Domain hits: {raw_summary['domain_hit_count']}")
    print(f"Unique domains: {raw_summary['unique_domains']}")
    print(
        f"Unique domains in at least {args.min_proteins_per_domain} proteins: {raw_summary['unique_domains_min_count']}"
    )
    print(f"Failures: {raw_summary['failure_count']}")

    print("\nHarmonized run results")
    print(f"Proteins: {harmonized_summary['protein_count']}")
    print(f"Domain hits: {harmonized_summary['domain_hit_count']}")
    print(f"Unique domains: {harmonized_summary['unique_domains']}")
    print(
        f"Unique domains in at least {args.min_proteins_per_domain} proteins: {harmonized_summary['unique_domains_min_count']}"
    )
    print(f"Failures: {harmonized_summary['failure_count']}")

    print("\nSource summary")
    for row in source_summary_rows:
        print(
            f"{row['source']}: hits={row['domain_hits']}, unique_raw={row['unique_raw_domain_keys']}, "
            f"mapped_hits={row['mapped_hits']}, unmapped_hits={row['unmapped_hits']}, ambiguous_hits={row['ambiguous_hits']}"
        )
    print(f"\nOutput directory: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
