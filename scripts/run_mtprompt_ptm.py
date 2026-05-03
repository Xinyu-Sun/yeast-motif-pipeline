#!/usr/bin/env python3
"""Prepare and normalize optional MTPrompt-PTM predictions.

This script intentionally does not vendor MTPrompt-PTM code, weights, or data.
It can generate per-type FASTA files from ``polymerase_tf_proteins.csv`` and
normalize prediction CSVs produced by an external MTPrompt checkout or Docker
workflow.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROTEINS = REPO_ROOT / "data" / "output" / "polymerase_tf_proteins.csv"
DEFAULT_KNOWN_PTMS = REPO_ROOT / "data" / "output" / "polymerase_tf_ptm_sites.csv"
DEFAULT_OUTDIR = REPO_ROOT / "data" / "output" / "mtprompt_ptm"
DEFAULT_PREDICTIONS_OUT = REPO_ROOT / "data" / "output" / "polymerase_tf_ptm_predictions_mtprompt.csv"
DEFAULT_CANDIDATES_OUT = REPO_ROOT / "data" / "output" / "polymerase_tf_ptm_candidate_sites.csv"

MTPROMPT_TYPES = [
    "Phosphorylation_S",
    "Phosphorylation_T",
    "Phosphorylation_Y",
    "Ubiquitination_K",
    "Acetylation_K",
    "OlinkedGlycosylation_S",
    "OlinkedGlycosylation_T",
    "Methylation_K",
    "Methylation_R",
    "NlinkedGlycosylation_N",
    "Palmitoylation_C",
    "Sumoylation_K",
    "Succinylation_K",
]

TYPE_TO_FAMILY = {
    "Phosphorylation_S": "phosphorylation",
    "Phosphorylation_T": "phosphorylation",
    "Phosphorylation_Y": "phosphorylation",
    "Ubiquitination_K": "ubiquitination",
    "Acetylation_K": "acetylation",
    "OlinkedGlycosylation_S": "other",
    "OlinkedGlycosylation_T": "other",
    "Methylation_K": "methylation",
    "Methylation_R": "methylation",
    "NlinkedGlycosylation_N": "other",
    "Palmitoylation_C": "palmitoylation",
    "Sumoylation_K": "sumoylation",
    "Succinylation_K": "succinylation",
}


def clean_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).replace("\xa0", " ").strip()


def normalized_header(value: str) -> str:
    return "".join(ch for ch in clean_text(value).lower() if ch.isalnum())


def int_or_none(value: object) -> int | None:
    try:
        if value == "":
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def float_or_blank(value: object) -> float | str:
    try:
        if value == "":
            return ""
        return float(value)
    except (TypeError, ValueError):
        return clean_text(value)


def read_table(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def protein_lookup(proteins: List[dict]) -> Dict[str, dict]:
    lookup: Dict[str, dict] = {}
    for row in proteins:
        aliases = {
            row.get("y_name", ""),
            row.get("standard_name", ""),
            row.get("protein", ""),
            row.get("systematic_name", ""),
            row.get("uniprot_id", ""),
            row.get("accession", ""),
        }
        header = row.get("fasta_header", "")
        if header:
            token = header.lstrip(">").split()[0]
            aliases.add(token)
            parts = token.split("|")
            if len(parts) >= 3:
                aliases.add(parts[2])
        for alias in aliases:
            if clean_text(alias):
                lookup[clean_text(alias)] = row
    return lookup


def write_fasta_inputs(proteins: List[dict], outdir: Path) -> List[Path]:
    fasta_dir = outdir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for ptm_type in MTPROMPT_TYPES:
        path = fasta_dir / f"{ptm_type}.fasta"
        with path.open("w", encoding="utf-8") as handle:
            for row in proteins:
                sequence = clean_text(row.get("amino_acid_sequence", "")).replace(" ", "").replace("\n", "")
                if not sequence:
                    continue
                header = row.get("y_name") or row.get("standard_name") or row.get("protein")
                handle.write(f">{header}\n")
                for idx in range(0, len(sequence), 80):
                    handle.write(sequence[idx : idx + 80] + "\n")
        written.append(path)
    return written


def run_external_commands(args: argparse.Namespace, fasta_paths: List[Path]) -> None:
    if args.external_command:
        for fasta_path in fasta_paths:
            command = [
                part.format(
                    fasta=str(fasta_path),
                    type=fasta_path.stem,
                    outdir=str(args.prediction_dir),
                    mtprompt_dir=str(args.mtprompt_dir or ""),
                )
                for part in args.external_command
            ]
            subprocess.run(command, check=True, cwd=args.mtprompt_dir or None)
    elif args.docker:
        args.prediction_dir.mkdir(parents=True, exist_ok=True)
        mounted_fasta = (args.outdir / "fasta").resolve()
        for fasta_path in fasta_paths:
            data_path = f"/app/data/{fasta_path.name}"
            subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--name",
                    f"mtprompt-{fasta_path.stem.lower()}",
                    "-v",
                    f"{mounted_fasta}:/app/data:ro",
                    "-v",
                    f"{args.prediction_dir.resolve()}:/app/predictions",
                    "hanye0311/mtprompt:v1",
                    "python3",
                    "test.py",
                    "--data_path",
                    data_path,
                    "--PTM_type",
                    fasta_path.stem,
                    "--save_path",
                    "/app/predictions",
                ],
                check=True,
            )


def first_value(row: dict, aliases: List[str]) -> str:
    normalized = {normalized_header(key): value for key, value in row.items()}
    for alias in aliases:
        value = clean_text(normalized.get(normalized_header(alias), ""))
        if value:
            return value
    return ""


def infer_type(path: Path, row: dict) -> str:
    explicit = first_value(row, ["ptm_type", "type", "modification", "modification_type", "task"])
    if explicit in MTPROMPT_TYPES:
        return explicit
    for ptm_type in MTPROMPT_TYPES:
        if ptm_type.lower() in path.stem.lower():
            return ptm_type
    return explicit


def normalize_prediction_rows(prediction_dir: Path, proteins: List[dict]) -> List[dict]:
    lookup = protein_lookup(proteins)
    rows: List[dict] = []
    for path in sorted(prediction_dir.glob("*.csv")):
        for raw in read_table(path):
            protein_key = first_value(raw, ["y_name", "protein", "protein_id", "id", "name", "sequence_id"])
            protein = lookup.get(protein_key)
            if protein is None:
                continue
            site = int_or_none(first_value(raw, ["site_index", "position", "pos", "site", "index"]))
            ptm_type = infer_type(path, raw)
            expected_residue = ptm_type.rsplit("_", 1)[-1] if "_" in ptm_type else first_value(raw, ["residue", "aa"])
            residue = (first_value(raw, ["site_residue", "residue", "aa", "amino_acid"]) or expected_residue).upper()
            sequence = clean_text(protein.get("amino_acid_sequence", "")).replace(" ", "").replace("\n", "")
            if site is None or site < 1 or site > len(sequence):
                continue
            observed_residue = sequence[site - 1].upper()
            if residue and observed_residue != residue:
                continue
            score = float_or_blank(first_value(raw, ["score", "probability", "prob", "prediction_score", "confidence"]))
            rows.append(
                {
                    "category": protein.get("category", ""),
                    "protein": protein.get("protein", ""),
                    "y_name": protein.get("y_name", ""),
                    "standard_name": protein.get("standard_name", ""),
                    "uniprot_id": protein.get("uniprot_id", ""),
                    "site_index": site,
                    "site_residue": observed_residue,
                    "mtp_type": ptm_type,
                    "ptm_family": TYPE_TO_FAMILY.get(ptm_type, "other"),
                    "score": score,
                    "source_file": str(path),
                }
            )
    rows.sort(key=lambda row: (row["category"], row["protein"], row["site_index"], row["mtp_type"]))
    return rows


def build_candidate_sites(predictions: List[dict], known_ptms: List[dict]) -> List[dict]:
    known_by_key: Dict[tuple[str, str, str, str], dict] = {}
    for row in known_ptms:
        key = (row.get("y_name", ""), clean_text(row.get("site_index")), row.get("site_residue", ""), row.get("ptm_family", ""))
        known_by_key[key] = row

    grouped: Dict[tuple[str, str, str, str], dict] = {}
    for key, row in known_by_key.items():
        grouped[key] = {
            "category": row.get("category", ""),
            "protein": row.get("protein", ""),
            "y_name": row.get("y_name", ""),
            "standard_name": row.get("standard_name", ""),
            "site_index": row.get("site_index", ""),
            "site_residue": row.get("site_residue", ""),
            "ptm_family": row.get("ptm_family", ""),
            "evidence_tier": "known_sgd",
            "known_sgd_ptm_evidence_ids": row.get("ptm_evidence_ids", ""),
            "mtp_types": set(),
            "prediction_count": 0,
            "max_score": "",
        }

    for row in predictions:
        key = (row["y_name"], clean_text(row["site_index"]), row["site_residue"], row["ptm_family"])
        group = grouped.setdefault(
            key,
            {
                "category": row.get("category", ""),
                "protein": row.get("protein", ""),
                "y_name": row.get("y_name", ""),
                "standard_name": row.get("standard_name", ""),
                "site_index": row.get("site_index", ""),
                "site_residue": row.get("site_residue", ""),
                "ptm_family": row.get("ptm_family", ""),
                "evidence_tier": "predicted_known_overlap" if key in known_by_key else "predicted_novel",
                "known_sgd_ptm_evidence_ids": known_by_key.get(key, {}).get("ptm_evidence_ids", ""),
                "mtp_types": set(),
                "prediction_count": 0,
                "max_score": "",
            },
        )
        if group["evidence_tier"] == "known_sgd":
            group["evidence_tier"] = "predicted_known_overlap"
        group["prediction_count"] += 1
        group["mtp_types"].add(row.get("mtp_type", ""))
        score = row.get("score", "")
        if isinstance(score, float) and (group["max_score"] == "" or score > group["max_score"]):
            group["max_score"] = score

    candidates: List[dict] = []
    for group in grouped.values():
        candidates.append(
            {
                **{k: v for k, v in group.items() if not isinstance(v, set)},
                "mtp_types": "; ".join(sorted(v for v in group["mtp_types"] if v)),
            }
        )
    candidates.sort(key=lambda row: (row["category"], row["protein"], int_or_none(row["site_index"]) or 0, row["ptm_family"]))
    return candidates


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, default=DEFAULT_PROTEINS)
    parser.add_argument("--known-ptms", type=Path, default=DEFAULT_KNOWN_PTMS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--prediction-dir", type=Path, default=DEFAULT_OUTDIR / "predictions")
    parser.add_argument("--predictions-out", type=Path, default=DEFAULT_PREDICTIONS_OUT)
    parser.add_argument("--candidates-out", type=Path, default=DEFAULT_CANDIDATES_OUT)
    parser.add_argument("--prepare-only", action="store_true", help="Only write MTPrompt FASTA inputs")
    parser.add_argument("--docker", action="store_true", help="Run Docker image hanye0311/mtprompt:v1 for each FASTA input")
    parser.add_argument("--mtprompt-dir", type=Path, help="External MTPrompt-PTM checkout directory")
    parser.add_argument(
        "--external-command",
        nargs="+",
        help="Command template to run per FASTA; placeholders: {fasta}, {type}, {outdir}, {mtprompt_dir}",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    proteins = read_table(args.proteins)
    fasta_paths = write_fasta_inputs(proteins, args.outdir)
    print(f"Wrote {len(fasta_paths)} FASTA files under {args.outdir / 'fasta'}")

    if args.prepare_only:
        return 0

    args.prediction_dir.mkdir(parents=True, exist_ok=True)
    run_external_commands(args, fasta_paths)
    predictions = normalize_prediction_rows(args.prediction_dir, proteins)
    known_ptms = read_table(args.known_ptms) if args.known_ptms.exists() else []
    candidates = build_candidate_sites(predictions, known_ptms)

    prediction_fields = [
        "category",
        "protein",
        "y_name",
        "standard_name",
        "uniprot_id",
        "site_index",
        "site_residue",
        "mtp_type",
        "ptm_family",
        "score",
        "source_file",
    ]
    candidate_fields = [
        "category",
        "protein",
        "y_name",
        "standard_name",
        "site_index",
        "site_residue",
        "ptm_family",
        "evidence_tier",
        "known_sgd_ptm_evidence_ids",
        "prediction_count",
        "mtp_types",
        "max_score",
    ]
    write_csv(args.predictions_out, predictions, prediction_fields)
    write_csv(args.candidates_out, candidates, candidate_fields)
    print(f"Wrote MTPrompt predictions: {args.predictions_out}")
    print(f"Wrote PTM candidate sites: {args.candidates_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
