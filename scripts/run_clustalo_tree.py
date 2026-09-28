#!/usr/bin/env python3
"""Align a protein CSV with Clustal Omega and write the tree used to order figure columns.

The runner can use either a local ``clustalo`` binary or the EMBL-EBI Job
Dispatcher REST service that backs the Clustal Omega web form.

Both backends write the same three files for the same FASTA:

- the alignment (Clustal format with residue numbers)
- the "phylotree": EBI's neighbour-joining tree computed from the finished
  alignment without distance corrections (the tree the web form offers as
  "Phylogenetic Tree")
- a percent identity matrix computed from the finished alignment

The local backend runs Clustal Omega with EBI's settings, then computes the
neighbour-joining tree and identity matrix the way EBI's ClustalW 2.1 step
does. In testing this reproduced EBI's alignment and tree files
byte-for-byte. Clustal Omega results depend on the order of sequences in the
FASTA, which is written in protein-CSV row order.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import time
import uuid
from itertools import combinations
from pathlib import Path
from urllib import error, parse, request
from xml.etree import ElementTree


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROTEINS = REPO_ROOT / "data" / "output" / "polymerase_tf_proteins.csv"
DEFAULT_FASTA = REPO_ROOT / "data" / "input" / "clustalo-extended.fasta"
DEFAULT_ALIGNMENT = REPO_ROOT / "data" / "input" / "clustalo-extended.clu"
DEFAULT_TREE = REPO_ROOT / "data" / "input" / "clustalo-extended.phylotree"
DEFAULT_DISTMAT = REPO_ROOT / "data" / "input" / "clustalo-extended.distmat"
EBI_REST_BASE = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
EBI_FINISHED_STATUS = "FINISHED"
EBI_FAILURE_STATUSES = {"ERROR", "FAILURE", "NOT_FOUND"}
EBI_RESULT_CANDIDATES = {
    "alignment": ("aln-clustal_num", "aln-clustal", "aln-clustalw", "aln"),
    "guide_tree": ("tree", "guidetree"),
    "phylo_tree": ("phylotree",),
    "distmat": ("distmat", "pim", "percentid"),
}


def clean_text(value: object) -> str:
    return "" if value is None else str(value).strip()


def fasta_id(row: dict) -> str:
    y_name = clean_text(row.get("y_name"))
    if y_name:
        return y_name
    label = clean_text(row.get("standard_name") or row.get("protein"))
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_")
    return label


def normalize_sequence(sequence: str) -> str:
    return re.sub(r"[^A-Za-z*]", "", clean_text(sequence)).upper()


def read_proteins(path: Path, category: str | None = None) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if category:
        rows = [row for row in rows if row.get("category") == category]
    if not rows:
        raise SystemExit(f"No protein rows found in {path}")
    return rows


def fasta_records(rows: list[dict]) -> list[str]:
    seen_ids: set[str] = set()
    missing: list[str] = []
    records: list[str] = []
    for row in rows:
        seq_id = fasta_id(row)
        if not seq_id:
            missing.append(clean_text(row.get("protein")) or "<missing protein>")
            continue
        if seq_id in seen_ids:
            raise SystemExit(f"Duplicate FASTA id for Clustal Omega: {seq_id}")
        seen_ids.add(seq_id)

        sequence = normalize_sequence(row.get("amino_acid_sequence", ""))
        if not sequence:
            missing.append(seq_id)
            continue
        description = clean_text(row.get("standard_name") or row.get("protein"))
        subgroup = clean_text(row.get("pol_i_subgroup"))
        if subgroup:
            description = f"{description} {subgroup}".strip()
        records.append(f">{seq_id} {description}".rstrip())
        records.extend(sequence[idx : idx + 80] for idx in range(0, len(sequence), 80))

    if missing:
        raise SystemExit("Missing amino_acid_sequence for: " + ", ".join(missing))

    return records


def write_fasta(path: Path, rows: list[dict]) -> str:
    records = fasta_records(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    fasta_text = "\n".join(records) + "\n"
    path.write_text(fasta_text, encoding="utf-8")
    return fasta_text


def clustalo_command(args: argparse.Namespace) -> list[str]:
    """Local Clustal Omega call matching the EBI full-distance settings and output format."""
    clustalo = shutil.which(args.clustalo_bin) or args.clustalo_bin
    cmd = [
        clustalo,
        "-i",
        str(args.fasta),
        "-o",
        str(args.alignment_out),
        "--outfmt=clu",
        "--resno",
        "--full",
        "--output-order=tree-order",
        "--force",
    ]
    if args.local_tree == "guide":
        cmd.extend(["--guidetree-out", str(args.out)])
    return cmd


def read_clustal_alignment(path: Path) -> list[tuple[str, str]]:
    """Return (name, aligned sequence) pairs in file order from a Clustal-format alignment."""
    sequences: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("CLUSTAL") or line[0].isspace():
            continue
        parts = line.split()
        if len(parts) >= 2 and re.fullmatch(r"[A-Za-z*.\-]+", parts[1]):
            sequences[parts[0]] = sequences.get(parts[0], "") + parts[1]
    if len(sequences) < 3:
        raise SystemExit(f"Need at least 3 aligned sequences to build a tree: {path}")
    return list(sequences.items())


def identity_and_distance(alignment: list[tuple[str, str]]) -> tuple[list[list[float | None]], list[list[float]]]:
    """Percent identity and uncorrected p-distance over columns where neither sequence has a gap.

    Pairs that share no ungapped column get identity None (EBI prints -nan) and
    distance 0.0. This follows ClustalW, whose tree EBI returns, even though it
    treats such pairs as identical.
    """
    sequences = [sequence for _, sequence in alignment]
    count = len(sequences)
    identity: list[list[float | None]] = [[100.0 if i == j else None for j in range(count)] for i in range(count)]
    distance = [[0.0] * count for _ in range(count)]
    for i, j in combinations(range(count), 2):
        compared = mismatches = 0
        for a, b in zip(sequences[i], sequences[j]):
            if a != "-" and b != "-":
                compared += 1
                mismatches += a != b
        if compared:
            identity[i][j] = identity[j][i] = 100.0 * (compared - mismatches) / compared
            distance[i][j] = distance[j][i] = mismatches / compared
    return identity, distance


def neighbor_joining_newick(names: list[str], distance: list[list[float]]) -> str:
    """Neighbour-joining tree in the same form ClustalW writes it.

    Uses Saitou and Nei's S criterion (not the algebraically equivalent Q form,
    which rounds differently at near-ties), keeps the first minimum in ClustalW's
    scan order, roots at the last three clusters, and orders children by their
    first sequence in alignment order.
    """
    count = len(names)
    d = [row[:] for row in distance]
    clusters: dict[int, tuple] = {i: (i,) for i in range(count)}
    alive = list(range(count))
    while len(alive) > 3:
        m = len(alive)
        row_sum = {i: sum(d[i][k] for k in alive if k != i) for i in alive}
        total = sum(d[i][k] for i, k in combinations(alive, 2))
        best = None
        for jj in alive:
            for ii in alive:
                if ii >= jj:
                    break
                dij = d[ii][jj]
                score = (
                    (row_sum[ii] + row_sum[jj] - 2 * dij) / (2 * (m - 2))
                    + dij / 2
                    + (total - row_sum[ii] - row_sum[jj] + dij) / (m - 2)
                )
                if best is None or score < best[0]:
                    best = (score, ii, jj)
        _, i, j = best
        length_i = 0.5 * d[i][j] + (row_sum[i] - row_sum[j]) / (2 * (m - 2))
        length_j = d[i][j] - length_i
        clusters[i] = ((clusters[i], length_i), (clusters[j], length_j))
        for k in alive:
            if k not in (i, j):
                d[i][k] = d[k][i] = 0.5 * (d[i][k] + d[j][k] - d[i][j])
        alive.remove(j)

    a, b, c = alive
    length_a = 0.5 * (d[a][b] + d[a][c] - d[b][c])
    root = [(clusters[a], length_a), (clusters[b], d[a][b] - length_a), (clusters[c], d[a][c] - length_a)]

    def first_index(cluster: tuple) -> int:
        return cluster[0] if len(cluster) == 1 else min(first_index(child) for child, _ in cluster)

    def render(children: list[tuple]) -> str:
        parts = []
        for cluster, length in sorted(children, key=lambda item: first_index(item[0])):
            label = names[cluster[0]] if len(cluster) == 1 else render(list(cluster)) + "\n"
            parts.append(f"{label}:{length:.5f}")
        return "(\n" + ",\n".join(parts) + ")"

    return render(root) + ";\n"


def format_percent_identity_matrix(names: list[str], identity: list[list[float | None]]) -> str:
    width = max(len(name) for name in names) + 1
    lines = ["#", "#", "#  Percent Identity  Matrix - computed from the alignment (ClustalW 2.1 convention)", "#", "#", ""]
    for index, name in enumerate(names, start=1):
        values = "".join("    -nan" if value is None else f"{value:8.2f}" for value in identity[index - 1])
        lines.append(f"{index:6d}: {name:<{width}}{values}")
    return "\n".join(lines) + "\n"


def write_tree_and_identity(args: argparse.Namespace) -> None:
    alignment = read_clustal_alignment(args.alignment_out)
    names = [name for name, _ in alignment]
    identity, distance = identity_and_distance(alignment)
    args.distmat_out.parent.mkdir(parents=True, exist_ok=True)
    args.distmat_out.write_text(format_percent_identity_matrix(names, identity), encoding="utf-8")
    print(f"Wrote percent identity matrix: {args.distmat_out}")
    no_overlap = sum(1 for i, j in combinations(range(len(names)), 2) if identity[i][j] is None)
    if no_overlap:
        print(f"Note: {no_overlap} sequence pair(s) share no ungapped alignment column; the tree treats them as distance 0.")
    if args.local_tree == "phylotree":
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(neighbor_joining_newick(names, distance), encoding="utf-8")
        print(f"Wrote tree (phylotree): {args.out}")


def ebi_url(*parts: str) -> str:
    return "/".join([EBI_REST_BASE.rstrip("/"), *(parse.quote(part.strip("/")) for part in parts)])


def form_request(url: str, fields: dict[str, str]) -> request.Request:
    boundary = f"----codex-clustalo-{uuid.uuid4().hex}"
    body_parts: list[bytes] = []
    for name, value in fields.items():
        body_parts.append(f"--{boundary}\r\n".encode("utf-8"))
        body_parts.append(
            (
                f'Content-Disposition: form-data; name="{name}"\r\n'
                "Content-Type: text/plain; charset=utf-8\r\n\r\n"
            ).encode("utf-8")
        )
        body_parts.append(value.encode("utf-8"))
        body_parts.append(b"\r\n")
    body_parts.append(f"--{boundary}--\r\n".encode("utf-8"))
    body = b"".join(body_parts)
    return request.Request(
        url,
        data=body,
        headers={
            "Content-Type": f"multipart/form-data; boundary={boundary}",
            "User-Agent": "yeast-motif-pipeline/1.0",
        },
        method="POST",
    )


def read_url(url: str, timeout: int) -> str:
    try:
        with request.urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"EBI request failed with HTTP {exc.code}: {detail}") from exc
    except error.URLError as exc:
        raise SystemExit(f"EBI request failed: {exc.reason}") from exc


def submit_ebi_job(args: argparse.Namespace, fasta_text: str) -> str:
    email = clean_text(args.email or os.environ.get("EBI_EMAIL"))
    if not email:
        raise SystemExit(
            "The EBI Clustal Omega REST service requires an email address. "
            "Pass --email or set EBI_EMAIL."
        )

    if args.ebi_mode == "web-default":
        mbed = "true"
        mbediteration = "true"
        dismatout = "false"
    else:
        mbed = "false"
        mbediteration = "false"
        dismatout = "true"

    fields = {
        "email": email,
        "title": clean_text(args.title) or args.fasta.stem,
        "stype": "protein",
        "sequence": fasta_text,
        "outfmt": args.ebi_outfmt,
        "order": args.ebi_order,
        "guidetreeout": "true",
        "dismatout": dismatout,
        "mbed": mbed,
        "mbediteration": mbediteration,
        "iterations": str(args.ebi_iterations),
        "gtiterations": str(args.ebi_gtiterations),
        "hmmiterations": str(args.ebi_hmmiterations),
    }

    req = form_request(ebi_url("run"), fields)
    try:
        with request.urlopen(req, timeout=args.timeout) as response:
            body = response.read().decode("utf-8").strip()
    except error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"EBI job submission failed with HTTP {exc.code}: {detail}") from exc
    except error.URLError as exc:
        raise SystemExit(f"EBI job submission failed: {exc.reason}") from exc

    if body.startswith("<?xml") or body.startswith("<"):
        raise SystemExit(f"EBI job submission failed: {body}")
    return body


def poll_ebi_job(job_id: str, poll_seconds: float, timeout: int, max_wait: float = 3600.0) -> str:
    status_url = ebi_url("status", job_id)
    deadline = time.monotonic() + max_wait
    while True:
        status = read_url(status_url, timeout).strip()
        print(f"EBI status: {status}")
        if status == EBI_FINISHED_STATUS:
            return status
        if status in EBI_FAILURE_STATUSES:
            raise SystemExit(f"EBI job ended with status: {status}")
        if time.monotonic() >= deadline:
            raise SystemExit(
                f"EBI job {job_id} still {status} after {max_wait:.0f} s. Check it later at {status_url} "
                "or rerun with a larger --max-wait."
            )
        time.sleep(poll_seconds)


def ebi_result_types(job_id: str, timeout: int) -> dict[str, str]:
    body = read_url(ebi_url("resulttypes", job_id), timeout)
    try:
        root = ElementTree.fromstring(body)
    except ElementTree.ParseError as exc:
        raise SystemExit(f"Could not parse EBI result types XML: {body}") from exc

    identifiers: dict[str, str] = {}
    for result_type in root.findall(".//type"):
        identifier = clean_text(result_type.findtext("identifier"))
        label = clean_text(result_type.findtext("label"))
        if identifier:
            identifiers[identifier] = label
    if not identifiers:
        raise SystemExit(f"EBI returned no result types: {body}")
    return identifiers


def choose_result_identifier(result_types: dict[str, str], kind: str) -> str | None:
    available = set(result_types)
    for candidate in EBI_RESULT_CANDIDATES[kind]:
        if candidate in available:
            return candidate
    for identifier in result_types:
        lowered = identifier.lower()
        if any(candidate in lowered for candidate in EBI_RESULT_CANDIDATES[kind]):
            return identifier
    return None


def download_ebi_result(job_id: str, result_id: str, out_path: Path, timeout: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    body = read_url(ebi_url("result", job_id, result_id), timeout)
    out_path.write_text(body, encoding="utf-8")


def run_ebi(args: argparse.Namespace, fasta_text: str) -> None:
    job_id = submit_ebi_job(args, fasta_text)
    print(f"EBI job id: {job_id}")
    if args.job_id_out:
        args.job_id_out.parent.mkdir(parents=True, exist_ok=True)
        args.job_id_out.write_text(job_id + "\n", encoding="utf-8")

    poll_ebi_job(job_id, args.poll_seconds, args.timeout, args.max_wait)
    result_types = ebi_result_types(job_id, args.timeout)
    print("EBI result types:")
    for identifier, label in sorted(result_types.items()):
        print(f"  {identifier}: {label}")

    alignment_id = choose_result_identifier(result_types, "alignment")
    tree_kind = "guide_tree" if args.ebi_tree_result == "guide" else "phylo_tree"
    tree_id = choose_result_identifier(result_types, tree_kind)
    distmat_id = choose_result_identifier(result_types, "distmat")

    if not alignment_id:
        raise SystemExit("EBI results did not include an alignment output")
    if not tree_id:
        raise SystemExit(f"EBI results did not include a {args.ebi_tree_result} tree output")

    download_ebi_result(job_id, alignment_id, args.alignment_out, args.timeout)
    download_ebi_result(job_id, tree_id, args.out, args.timeout)
    print(f"Wrote alignment: {args.alignment_out}")
    print(f"Wrote tree ({args.ebi_tree_result}): {args.out}")

    if args.ebi_mode == "full-distance":
        if not distmat_id:
            raise SystemExit("EBI results did not include a distance matrix output")
        download_ebi_result(job_id, distmat_id, args.distmat_out, args.timeout)
        print(f"Wrote distance matrix: {args.distmat_out}")
    elif distmat_id:
        download_ebi_result(job_id, distmat_id, args.distmat_out, args.timeout)
        print(f"Wrote distance matrix: {args.distmat_out}")
    else:
        print("No distance matrix requested for EBI web-default mode")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, default=DEFAULT_PROTEINS)
    parser.add_argument("--category", choices=["Pol I", "Pol II", "Pol III"], help="Optional broad category filter")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--alignment-out", type=Path, default=DEFAULT_ALIGNMENT)
    parser.add_argument("--out", type=Path, default=DEFAULT_TREE)
    parser.add_argument("--distmat-out", type=Path, default=DEFAULT_DISTMAT)
    parser.add_argument("--backend", choices=["local", "ebi"], default="local")
    parser.add_argument("--clustalo-bin", default="clustalo")
    parser.add_argument(
        "--local-tree",
        choices=["phylotree", "guide"],
        default="phylotree",
        help=(
            "Tree written to --out by the local backend. phylotree (default) is the same "
            "neighbour-joining tree EBI returns; guide is Clustal Omega's own guide tree."
        ),
    )
    parser.add_argument("--write-fasta-only", action="store_true", help="Write FASTA and print the Clustal command")
    parser.add_argument("--email", help="Email address required by the EBI Job Dispatcher REST service")
    parser.add_argument("--title", help="Optional EBI job title")
    parser.add_argument(
        "--ebi-mode",
        choices=["full-distance", "web-default"],
        default="full-distance",
        help=(
            "EBI parameter bundle. full-distance disables mBed and requests a distance "
            "matrix; web-default keeps the web form defaults and does not request one."
        ),
    )
    parser.add_argument("--ebi-outfmt", default="clustal_num", help="EBI alignment output format")
    parser.add_argument("--ebi-order", choices=["aligned", "input"], default="aligned")
    parser.add_argument(
        "--ebi-tree-result",
        choices=["phylotree", "guide"],
        default="phylotree",
        help="EBI tree result to save to --out. Use phylotree for figure-style trees; guide saves the Clustal guide tree.",
    )
    parser.add_argument("--ebi-iterations", type=int, default=0)
    parser.add_argument("--ebi-gtiterations", type=int, default=-1)
    parser.add_argument("--ebi-hmmiterations", type=int, default=-1)
    parser.add_argument("--poll-seconds", type=float, default=5.0)
    parser.add_argument("--timeout", type=int, default=60, help="Seconds to wait for each HTTP request")
    parser.add_argument("--max-wait", type=float, default=3600.0, help="Seconds to wait for the EBI job to finish")
    parser.add_argument("--job-id-out", type=Path, help="Optional file for the EBI job id")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = read_proteins(args.proteins, category=args.category)
    fasta_text = write_fasta(args.fasta, rows)
    cmd = clustalo_command(args)

    print(f"Wrote FASTA: {args.fasta}")
    print(f"Proteins: {len(rows)}")
    if args.backend == "local":
        print("Clustal Omega command:")
        print(" ".join(cmd))
    else:
        print("Clustal Omega backend: EMBL-EBI Job Dispatcher REST")
        print(f"EBI mode: {args.ebi_mode}")

    if args.write_fasta_only:
        return 0

    if args.backend == "ebi":
        run_ebi(args, fasta_text)
        return 0

    if not shutil.which(args.clustalo_bin) and not Path(args.clustalo_bin).exists():
        raise SystemExit(
            "Clustal Omega binary not found. Install `clustalo` or pass --clustalo-bin, "
            "then rerun this command."
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.alignment_out.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    if result.returncode != 0:
        raise SystemExit(result.returncode)

    print(f"Wrote alignment: {args.alignment_out}")
    if args.local_tree == "guide":
        print(f"Wrote tree (Clustal Omega guide tree): {args.out}")
    write_tree_and_identity(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
