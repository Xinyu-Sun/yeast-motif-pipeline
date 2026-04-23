# Yeast Motif Pipeline

A clean, GitHub-ready pipeline for extracting yeast polymerase transcription factor motif/domain annotations (SGD + InterPro), summarizing shared domains, and producing Clustal-Omega-ordered motif dot-matrix figures.

This repository adapts the original scripts into a reusable structure with relative paths and a safer plotting interface (CLI + JSON style config) so layout changes do not require code edits.

## Repository Layout

```text
yeast-motif-pipeline/
├── configs/
│   ├── plot_layout.example.json
│   └── plot_layout.compact.example.json
├── data/
│   ├── input/
│   │   ├── README.md
│   │   └── clustalo-all48.phylotree
│   └── output/
│       └── README.md
├── docs/
├── scripts/
│   ├── extract_sgd_domains.py
│   ├── summarize_shared_domains.py
│   └── plot_clustalo_motif_figure.py
├── .gitignore
└── README.md
```

## Requirements

- Python 3.9+
- Internet access for `extract_sgd_domains.py` (calls SGD + InterPro APIs)

No non-stdlib Python packages are required.

## Inputs

### 1) Protein list input (`.csv`, preferred)

The extractor can read a simple CSV table of proteins to look up:

```csv
category,protein,y_name,accession,fasta_header,amino_acid_sequence,raw_fasta
Pol II,Example TF 1,YAL001C,,,,
```

Required lookup column:
- Prefer `y_name`: SGD locus/systematic name used for API lookup.
- A single `protein`, `name`, or `standard_name` column also works if its values are valid SGD lookup names.

Optional columns:
- `category`: group label such as `Pol I`, `Pol II`, or `Pol III`
- `protein`: display name
- `accession`: accession identifier
- `fasta_header`, `amino_acid_sequence`, or `raw_fasta`: optional sequence metadata

Accepted aliases include `Y Name`, `systematic_name`, or `locus` for `y_name`;
`name` or `standard_name` for `protein`; and `sequence` or `aa_sequence` for
`amino_acid_sequence`.

Example template:
- `data/input/proteins.example.csv`

Run with `--csv`:

```bash
python3 scripts/extract_sgd_domains.py \
  --csv data/input/proteins.csv \
  --outdir data/output
```

### 2) Motif collection workbook (`.xlsx`, legacy)

The extractor expects the same column block layout as the original workbook:
- Pol I block: columns `A:D`
- Pol II block: columns `F:I`
- Pol III block: columns `K:N`
- Combined block: columns `P:S`

Default expected location:
- `data/input/polymerase_tf_sequences.xlsx`

Use `--xlsx` to pass a different file.

### 3) Clustal Omega tree (`.phylotree`)

Default tree:
- `data/input/clustalo-all48.phylotree`

This is the Newick guide tree used for protein ordering and dendrogram rendering.

## End-to-End Usage

Run all commands from repository root.

### Step 1: Extract domains from SGD + InterPro harmonization

```bash
python3 scripts/extract_sgd_domains.py \
  --csv data/input/proteins.csv \
  --outdir data/output
```

For the original workbook format, use `--xlsx data/input/polymerase_tf_sequences.xlsx`.

Useful options:
- `--force-refresh`: re-fetch API responses instead of using cache
- `--min-proteins-per-domain 2`: control `minN` matrix outputs
- `--sleep-seconds 0.15`: API pacing
- `--timeout 30`

Primary outputs include:
- `data/output/polymerase_tf_proteins.csv`
- `data/output/polymerase_tf_domain_hits_raw.csv`
- `data/output/polymerase_tf_domain_hits_harmonized.csv`
- `data/output/polymerase_tf_domain_matrix_*.csv`
- `data/output/run_summary.json`

### Step 2: Generate shared-domain category summary table

```bash
python3 scripts/summarize_shared_domains.py \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --out data/output/all48_shared_domain_category_summary.csv \
  --min-proteins 2
```

This table reports, per shared raw domain key:
- proteins sharing count
- category counts and percentages
- dominant category

### Step 3: Plot all proteins (default figure)

```bash
python3 scripts/plot_clustalo_motif_figure.py \
  --tree data/input/clustalo-all48.phylotree \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --proteins data/output/polymerase_tf_proteins.csv \
  --out data/output/all48_motif_figure.svg \
  --order-csv data/output/all48_motif_order.csv \
  --title "All 48 Protein TFs" \
  --subtitle "Hierarchical Clustering" \
  --min-proteins 2
```

Defaults:
- MobiDBLite is excluded
- legend is enabled

## Category-Specific Plotting

You can render Pol I / Pol II / Pol III subsets from the same combined tree (the script prunes the tree automatically).

### Example: Pol II only

```bash
python3 scripts/plot_clustalo_motif_figure.py \
  --category "Pol II" \
  --tree data/input/clustalo-all48.phylotree \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --proteins data/output/polymerase_tf_proteins.csv \
  --out data/output/pol2_motif_figure.svg \
  --order-csv data/output/pol2_motif_order.csv \
  --title "Pol II TFs" \
  --subtitle "Shared Domains (>=2 proteins)"
```

## Layout and Style Tweaking (No Code Edits)

Use either a JSON style config or direct CLI flags.

### Option A: JSON config (recommended for repeatability)

Base template:
- `configs/plot_layout.example.json`

Compact example:
- `configs/plot_layout.compact.example.json`

Run with config:

```bash
python3 scripts/plot_clustalo_motif_figure.py \
  --style-config configs/plot_layout.compact.example.json \
  --out data/output/all48_motif_figure_compact.svg
```

Generate a fresh default template:

```bash
python3 scripts/plot_clustalo_motif_figure.py \
  --write-default-style-config configs/my_style.json
```

### Option B: Direct CLI overrides

```bash
python3 scripts/plot_clustalo_motif_figure.py \
  --cell-width 12 \
  --cell-height 12 \
  --margin-left 180 \
  --margin-right 24 \
  --title-font-size 20 \
  --subtitle-font-size 13 \
  --title-subtitle-gap 16 \
  --subtitle-matrix-gap 20 \
  --legend-position top-right \
  --dot-fill "#2e7d32" \
  --dot-stroke "#1b5e20" \
  --category-color "Pol I=#0b69a3" \
  --category-color "Pol II=#ad2139" \
  --category-color "Pol III=#2d7c2f" \
  --out data/output/all48_motif_figure_tuned.svg
```

### Plotting controls available

- `--min-proteins`: shared-domain threshold
- `--exclude-mobidblite` / `--include-mobidblite`
- Margins: `--margin-left/right/top/bottom`
- Cell geometry: `--cell-width`, `--cell-height`
- Vertical spacing: `--title-subtitle-gap`, `--subtitle-matrix-gap`, `--dendrogram-gap`
- Typography: `--font-family`, `--title-font-size`, `--subtitle-font-size`, `--protein-label-font-size`, `--domain-label-font-size`
- Colors: `--dot-fill`, `--dot-stroke`, `--grid-color`, `--background-color`, `--category-color`
- Legend: `--no-legend`, `--legend-position`, `--legend-x`, `--legend-y`

## Output and Caching Notes

- All generated tables/figures are intended to live in `data/output/`.
- API responses are cached under `data/output/cache/`.
- `.gitignore` excludes generated outputs and common large assets from version control.

## Caveats

- SGD/InterPro responses can change over time, so exact row counts may drift.
- If your workbook schema differs from the expected columns, extraction will fail until the mapping is updated.
- SVG output is native from the plotting script. PNG/PDF export can be done separately if needed, but large generated images/PDFs should not be committed.
