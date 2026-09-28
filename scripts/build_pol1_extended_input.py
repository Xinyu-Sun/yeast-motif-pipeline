#!/usr/bin/env python3
"""Build an all-48-plus-Pol-I axial/peri-axial protein CSV.

The base 48-protein CSV already contains sequences. The Pol I peri-axial table adds
Pol I-associated peri-axial proteins but does not include FASTA sequences,
so this helper resolves those rows through UniProt and writes a combined CSV
that can be passed to ``extract_sgd_domains.py`` and the Clustal Omega runner.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from extract_sgd_domains import XlsxSheet, clean_text


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BASE_CSV = REPO_ROOT / "data" / "input" / "proteins.all48.generated.csv"
DEFAULT_OUT_CSV = REPO_ROOT / "data" / "input" / "proteins.all48_plus_pol1_axial_periaxial.csv"
DEFAULT_CACHE_DIR = REPO_ROOT / "data" / "output" / "cache" / "uniprot"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

FIELDNAMES = [
    "category",
    "protein",
    "y_name",
    "accession",
    "pol_i_subgroup",
    "source_table",
    "source_row",
    "source_note",
    "fasta_header",
    "amino_acid_sequence",
    "raw_fasta",
]


def normalize_gene_name(protein: str) -> str:
    first = clean_text(protein).split("/", 1)[0]
    return re.sub(r"[^A-Za-z0-9-]", "", first).upper()


def cache_name(gene: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", gene.upper()) + ".json"


def collect_gene_names(record: dict) -> tuple[str, str, list[str]]:
    standard_name = ""
    systematic_name = ""
    names: list[str] = []
    for gene in record.get("genes", []) or []:
        gene_name = clean_text((gene.get("geneName") or {}).get("value"))
        if gene_name and not standard_name:
            standard_name = gene_name
        if gene_name:
            names.append(gene_name)
        for item in gene.get("orderedLocusNames", []) or []:
            value = clean_text(item.get("value"))
            if value and not systematic_name:
                systematic_name = value
            if value:
                names.append(value)
        for item in gene.get("synonyms", []) or []:
            value = clean_text(item.get("value"))
            if value:
                names.append(value)
    return standard_name, systematic_name, names


def uniprot_query(gene: str) -> str:
    query = f"gene_exact:{gene} AND organism_id:559292 AND reviewed:true"
    fields = "accession,id,gene_names,sequence"
    return f"{UNIPROT_SEARCH}?query={quote(query)}&fields={quote(fields)}&format=json&size=5"


def fetch_uniprot_records(gene: str, cache_dir: Path, offline: bool, timeout: float) -> list[dict]:
    cache_path = cache_dir / cache_name(gene)
    if cache_path.exists():
        data = json.loads(cache_path.read_text(encoding="utf-8"))
    else:
        if offline:
            raise SystemExit(f"Missing cached UniProt lookup for {gene}: {cache_path}")
        request = Request(uniprot_query(gene), headers={"User-Agent": "yeast-motif-pipeline/1.0"})
        try:
            with urlopen(request, timeout=timeout) as response:
                data = json.loads(response.read().decode("utf-8"))
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
            raise SystemExit(f"Could not fetch UniProt record for {gene}: {exc}") from exc
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")

    results = data.get("results", []) if isinstance(data, dict) else []
    if not results:
        raise SystemExit(f"No reviewed UniProt result found for {gene}")
    return results


def choose_uniprot_hit(query_gene: str, results: list[dict]) -> tuple[dict, str]:
    """Pick the hit that really is ``query_gene``.

    UniProt's gene_exact search also matches aliases, so a name can return a
    different gene that lists it as a synonym. Prefer a hit whose main gene
    name or locus is the query. Accept an alias-only match only when exactly one
    hit has it (and say so in the returned note). Otherwise stop.
    """
    query = query_gene.upper()
    primary, alias = [], []
    for record in results:
        standard_name, systematic_name, names = collect_gene_names(record)
        if query in {standard_name.upper(), systematic_name.upper()}:
            primary.append(record)
        elif query in {name.upper() for name in names}:
            alias.append((record, standard_name, systematic_name))
    if primary:
        return primary[0], ""
    if len(alias) == 1:
        record, standard_name, systematic_name = alias[0]
        return record, f"Workbook name {query_gene} matched UniProt alias of {standard_name} ({systematic_name})"
    found = ", ".join(collect_gene_names(record)[0] or "?" for record in results)
    reason = "several genes list it as an alias" if alias else "no hit lists it as a gene name"
    raise SystemExit(f"Cannot resolve {query_gene} in UniProt: {reason} (hits: {found}). Check the workbook name.")


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[idx : idx + width] for idx in range(0, len(sequence), width))


def fasta_header(record: dict, standard_name: str) -> str:
    accession = clean_text(record.get("primaryAccession"))
    entry_id = clean_text(record.get("uniProtkbId")) or f"{standard_name.upper()}_YEAST"
    return f">sp|{accession}|{entry_id}"


def read_base_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"No rows found in base CSV: {path}")

    out: list[dict] = []
    for row in rows:
        normalized = {field: row.get(field, "") for field in FIELDNAMES}
        normalized.update(
            {
                "category": row.get("category", ""),
                "protein": row.get("protein", ""),
                "y_name": row.get("y_name", ""),
                "accession": row.get("accession", ""),
                "fasta_header": row.get("fasta_header", ""),
                "amino_acid_sequence": row.get("amino_acid_sequence", ""),
                "raw_fasta": row.get("raw_fasta", ""),
                "source_table": row.get("source_table", "") or "all48",
            }
        )
        if normalized["category"] == "Pol I" and not normalized["pol_i_subgroup"]:
            normalized["pol_i_subgroup"] = "axial"
        out.append(normalized)
    return out


def iter_pol_i_workbook_rows(path: Path) -> list[dict]:
    sheet = XlsxSheet(path)
    rows: list[dict] = []
    for row_num in range(3, sheet.max_row + 1):
        protein = sheet.get("A", row_num)
        y_name = sheet.get("B", row_num)
        if not protein:
            continue
        rows.append(
            {
                "source_row": str(row_num),
                "protein": protein,
                "workbook_y_name": y_name,
                "length_aa": sheet.get("C", row_num),
            }
        )
    if not rows:
        raise SystemExit(f"No Pol I protein rows found in workbook: {path}")
    return rows


def build_pol_i_row(source: dict, cache_dir: Path, offline: bool, timeout: float) -> dict:
    query_gene = normalize_gene_name(source["protein"])
    results = fetch_uniprot_records(query_gene, cache_dir=cache_dir, offline=offline, timeout=timeout)
    record, match_note = choose_uniprot_hit(query_gene, results)
    standard_name, systematic_name, _names = collect_gene_names(record)
    sequence = clean_text((record.get("sequence") or {}).get("value"))
    accession = clean_text(record.get("primaryAccession"))
    if not standard_name:
        standard_name = query_gene
    if not systematic_name:
        systematic_name = clean_text(source.get("workbook_y_name"))
    if not sequence or not accession or not systematic_name:
        raise SystemExit(f"UniProt record for {query_gene} is missing accession, locus, or sequence")

    source_note = match_note
    workbook_y_name = clean_text(source.get("workbook_y_name"))
    if workbook_y_name and workbook_y_name != systematic_name:
        correction = f"Workbook Y name {workbook_y_name} corrected to {systematic_name} by UniProt lookup for {standard_name}"
        source_note = f"{source_note}; {correction}" if source_note else correction

    try:
        workbook_length = int(float(source.get("length_aa") or 0))
    except ValueError:
        workbook_length = 0
    if workbook_length and workbook_length != len(sequence):
        length_note = f"Workbook length {workbook_length} aa differs from UniProt sequence length {len(sequence)} aa"
        source_note = f"{source_note}; {length_note}" if source_note else length_note

    header = fasta_header(record, standard_name)
    return {
        "category": "Pol I",
        "protein": standard_name,
        "y_name": systematic_name,
        "accession": accession,
        "pol_i_subgroup": "periaxial",
        "source_table": "Pol I peri-axial table",
        "source_row": source["source_row"],
        "source_note": source_note,
        "fasta_header": header,
        "amino_acid_sequence": sequence,
        "raw_fasta": f"{header}\n{wrap_sequence(sequence)}",
    }


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDNAMES})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-csv", type=Path, default=DEFAULT_BASE_CSV)
    parser.add_argument("--pol-i-workbook", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_CSV)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--offline", action="store_true", help="Use only cached UniProt lookups")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--sleep-seconds", type=float, default=0.1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.base_csv.exists():
        raise SystemExit(f"Base CSV not found: {args.base_csv}")
    if not args.pol_i_workbook.exists():
        raise SystemExit(f"Pol I workbook not found: {args.pol_i_workbook}")

    rows = read_base_rows(args.base_csv)
    seen_y_names = {row["y_name"] for row in rows if row.get("y_name")}
    additions: list[dict] = []
    workbook_rows = iter_pol_i_workbook_rows(args.pol_i_workbook)
    for idx, source_row in enumerate(workbook_rows, start=1):
        row = build_pol_i_row(source_row, cache_dir=args.cache_dir, offline=args.offline, timeout=args.timeout)
        if row["y_name"] in seen_y_names:
            print(f"Skipping duplicate Pol I row already present in base CSV: {row['protein']} ({row['y_name']})")
            continue
        seen_y_names.add(row["y_name"])
        additions.append(row)
        if args.sleep_seconds and idx < len(workbook_rows) and not args.offline:
            time.sleep(args.sleep_seconds)

    combined = rows + additions
    write_rows(args.out, combined)

    print(f"Wrote combined CSV: {args.out}")
    print(f"Base proteins: {len(rows)}")
    print(f"Added Pol I-associated proteins: {len(additions)}")
    print(f"Total proteins: {len(combined)}")
    print(f"Category counts: {dict(Counter(row['category'] for row in combined))}")
    print(
        "Pol I subgroup counts: "
        f"{dict(Counter(row['pol_i_subgroup'] for row in combined if row['category'] == 'Pol I'))}"
    )
    notes = [row for row in additions if row.get("source_note")]
    if notes:
        print("Source notes:")
        for row in notes:
            print(f"- {row['protein']} ({row['y_name']}): {row['source_note']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
