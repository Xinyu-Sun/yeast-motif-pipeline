# Yeast Motif Pipeline

A clean, GitHub-ready pipeline for extracting yeast polymerase transcription factor motif/domain annotations (SGD + InterPro), summarizing shared domains, and producing Clustal-Omega-ordered motif dot-matrix figures.

This repository adapts the original scripts into a reusable structure with relative paths and a safer plotting interface (CLI + JSON style config) so layout changes do not require code edits.

What changed in each update is listed in [`CHANGELOG.md`](CHANGELOG.md).

## Repository Layout

```text
yeast-motif-pipeline/
├── configs/
│   ├── plot_layout.example.json
│   └── plot_layout.compact.example.json
├── data/
│   ├── input/
│   │   ├── README.md
│   │   ├── clustalo-all48.phylotree
│   │   └── proteins.example.csv
│   └── output/
│       └── README.md
├── docs/
│   ├── clustalo_ebi_reproducibility.md
│   ├── matplotlib_renderer.md
│   └── pipeline_reference.md
├── scripts/
│   ├── build_output_index.py                    # HTML index of a run's figures and tables
│   ├── build_pol1_extended_input.py
│   ├── extract_sgd_domains.py
│   ├── plot_clustalo_motif_figure.py            # native SVG renderer (shared plotting code)
│   ├── plot_clustalo_motif_figure_matplotlib.py # recommended renderer
│   ├── plot_mobidblite_linear_schematic.py
│   ├── run_clustalo_tree.py
│   ├── run_mtprompt_ptm.py
│   ├── run_pipeline.py                          # one command: extract -> plots -> index.html
│   └── summarize_shared_domains.py
├── tests/
├── .gitignore
├── CHANGELOG.md
├── LICENSE
├── README.md
└── requirements-plotting.txt
```

## Requirements

- Python 3.9+
- Internet access for `extract_sgd_domains.py` (calls SGD + InterPro APIs)
- Clustal Omega (`clustalo`) to regenerate sequence-similarity guide trees
- Matplotlib for figure generation

No non-stdlib Python packages are required for extraction or summarization.

## Setup

Create a virtual environment and install the plotting dependency:

macOS / Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements-plotting.txt
```

Windows (PowerShell):

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements-plotting.txt
```

The plotting requirements currently install Matplotlib. Extraction and summary
scripts use only the Python standard library.

### Running the commands on Windows

The commands in this README are written for macOS/Linux shells. On Windows:

- Use `python` (or `py`) instead of `python3`.
- Lines ending in `\` are continued on the next line. In PowerShell, replace the
  trailing `\` with a backtick (`` ` ``); in Command Prompt, use `^`; or put the
  whole command on one line.
- Forward slashes in paths (`data/output/...`) work on Windows too.
- Quote paths that contain spaces, e.g. `--out "data/output/my run/figure.svg"`.
- If you edit input CSVs in Excel, save them as "CSV UTF-8". The scripts accept
  the byte-order mark that Excel adds.
- Native-renderer PDF/PNG export uses `rsvg-convert` when it is installed
  (common on macOS via Homebrew, rarely on Windows). The Matplotlib renderer
  does not need it, so it is the simplest choice on Windows.

## Quick Start: Whole Pipeline in One Command

`scripts/run_pipeline.py` runs every step below for one protein CSV and ends
by writing `index.html` in the output folder. The command is the same on macOS
and Windows (on Windows, use `python` instead of `python3`).

```bash
python3 scripts/run_pipeline.py --csv data/input/proteins.csv --tree data/input/my_tree.phylotree --outdir data/output/my_run
```

It writes, into `--outdir`:

- the extraction tables (Step 1) and the shared-domain summary (Step 2)
- dot matrices for all proteins and for each of Pol I / Pol II / Pol III that
  is present, named `<prefix>_<scope>_<raw|deduped>_by_<tree|motifs>`:
  - `raw`: motif IDs as reported by each source database
  - `deduped`: collapsed to InterPro entries, empty columns hidden
  - `by_tree`: columns in Clustal Omega tree order
  - `by_motifs`: columns clustered by shared motifs (see "Column order" below)

  Each figure has SVG, PNG and PDF files and a `_columns.csv` listing the
  column order.
- the MobiDBLite linear schematic (Step 4)
- `index.html` (Step 5)

Other options:

- `--prefix`: file-name prefix for the figures.
- `--order-by tree|motifs|both`: which column orders to draw; the default is
  `both`. With `--order-by motifs` no tree file is needed.
- `--skip-extract`: replot without calling SGD/InterPro again.
- `--renderer native`: use the native SVG renderer.

SGD and InterPro annotations change over time, so a later extraction can
differ slightly. Keep the output folder's `cache/` to reuse earlier responses.

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
- `pol_i_subgroup`: optional Pol I subgroup such as `axial` or `periaxial`
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

To build the 48-protein input plus the Pol I peri-axial table:

```bash
python3 scripts/build_pol1_extended_input.py \
  --base-csv data/input/proteins.all48.generated.csv \
  --pol-i-workbook "/path/to/pol_i_periaxial_table.xlsx" \
  --out data/input/proteins.pol1_extended.csv
```

The helper labels the original Pol I rows as `axial` and the newly added
workbook rows as `periaxial` for the plot label `Pol I peri-axial`. It also fills
missing FASTA sequences from UniProt and records any corrected workbook locus in
`source_note`. It uses a UniProt hit only if that hit really is the named gene
(see `choose_uniprot_hit` in the script).

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

This Newick tree sets the protein (column) order and the dendrogram in the
dot-matrix figures.

`scripts/run_clustalo_tree.py` makes a tree from any protein CSV that contains
`amino_acid_sequence`. It writes four files:

- `--fasta`: the sequences, in CSV row order
- `--alignment-out`: the Clustal Omega alignment
- `--out`: the tree
- `--distmat-out`: a percent identity matrix

The tree is the one EBI's Clustal Omega page offers as "Phylogenetic Tree": a
neighbour-joining tree computed from the finished alignment.

Local Clustal Omega and the EBI web service give identical results for the
same FASTA. In testing, the local run reproduced EBI's alignment and tree files
byte-for-byte (see `docs/clustalo_ebi_reproducibility.md`).
Choose whichever is easier:

- **Local** (needs `clustalo` installed; no internet or email):

```bash
python3 scripts/run_clustalo_tree.py \
  --proteins data/output/polymerase_tf_proteins.csv \
  --fasta data/input/clustalo-extended.fasta \
  --alignment-out data/input/clustalo-extended.clu \
  --out data/input/clustalo-extended.phylotree \
  --distmat-out data/input/clustalo-extended.distmat
```

- **EBI web service** (nothing to install; needs internet and an email address
  for the Job Dispatcher REST API):

```bash
python3 scripts/run_clustalo_tree.py \
  --backend ebi \
  --email name@example.edu \
  --proteins data/output/polymerase_tf_proteins.csv \
  --fasta data/input/clustalo-extended.fasta \
  --alignment-out data/input/clustalo-extended.clu \
  --out data/input/clustalo-extended.phylotree \
  --distmat-out data/input/clustalo-extended.distmat \
  --ebi-tree-result phylotree \
  --job-id-out data/input/clustalo-extended.ebi-job.txt
```

The default EBI mode is `--ebi-mode full-distance`, which disables mBed and
requests the identity matrix. `--ebi-mode web-default` keeps the web form's
mBed default. In testing, both modes gave the same alignment and tree.

For Clustal Omega's own guide tree instead of the neighbour-joining tree, use
`--local-tree guide` (local) or `--ebi-tree-result guide` (EBI).

**Two cautions when reading the tree:**

- **Sequence order in the FASTA changes the result.** Clustal Omega is not
  order-independent. The script writes the FASTA in CSV row order, so keep the
  row order fixed. A tree made by pasting the same sequences into the web form
  in a different order will not match. In testing, reshuffled orders gave trees
  barely closer to the original than random trees, even ignoring branch flips.
- **Transcription factors from different families are mostly unrelated by
  sequence.** In testing on a lab dataset, 96% of pairs shared under 20% identity.
  Pairs with no aligned residue in common are treated as identical by the
  tree step. For such sets, treat the dendrogram as a display order, not as
  evidence of relationships. Ordering by shared motifs (`--order-by motifs`,
  below) is often the more useful view.

## End-to-End Usage

Run all commands from repository root.

### Step 1: Extract domains from SGD + InterPro harmonization

Start from `data/input/proteins.example.csv` (copy it to `data/input/proteins.csv`
and fill in your proteins), or pass your own CSV path.

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
- `data/output/polymerase_tf_ptm_sites_raw.csv`
- `data/output/polymerase_tf_ptm_sites.csv`
- `data/output/polymerase_tf_ptm_domain_intersections.csv`
- `data/output/polymerase_tf_ptm_summary_by_protein.csv`
- `data/output/run_summary.json`

PTM options:
- `--skip-ptms`: disable SGD PTM fetching and PTM table outputs.
- `--ptm-boundary-window 10`: control the amino-acid window used to call a PTM near a domain boundary.

SGD known PTMs are treated as the authoritative PTM evidence layer. Known PTM
records are cached separately under `data/output/cache/sgd/*.ptms.json`.

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

### Step 3: Plot all proteins with Matplotlib

```bash
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
  --tree data/input/clustalo-all48.phylotree \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --proteins data/output/polymerase_tf_proteins.csv \
  --out data/output/all48_motif_figure.svg \
  --png data/output/all48_motif_figure.png \
  --pdf data/output/all48_motif_figure.pdf \
  --order-csv data/output/all48_motif_order.csv \
  --min-proteins 2
```

Defaults:
- MobiDBLite is excluded
- legend is enabled
- SVG is written by default; PNG/PDF are optional exports
- PTM tables are not rendered on motif/domain matrix plots
- title "All proteins" (or "Pol II proteins" with `--category "Pol II"`); the
  subtitle states the column order

#### Column order: tree or shared motifs

`--order-by` chooses how protein columns are ordered and which dendrogram is
drawn under the matrix. Each choice has its own default subtitle, and a
suggested figure legend is recorded in the figure's `.plot.json` file and
shown in `index.html`.

- `--order-by tree` (default): leaf order of the Clustal Omega tree given by
  `--tree`.
  - Subtitle: "Ordered by Clustal Omega sequence similarity (neighbour-joining)".
  - Suggested legend: "Columns follow a neighbour-joining tree built from a
    Clustal Omega alignment. When most proteins share under 20% sequence
    identity, the dendrogram is a display order, not evidence of evolutionary
    relationships."
- `--order-by motifs`: hierarchical clustering on the plotted dot matrix
  itself. No tree file is needed.
  - Two proteins' distance is 1 minus the fraction of their motifs they share
    (Jaccard), with average linkage.
  - Proteins with the same motifs end up side by side.
  - Proteins with no plotted motif are placed last.
  - The order does not depend on the order of rows in the input CSV.
  - Subtitle: "Ordered by shared motifs (Jaccard distance, average linkage)".

### Step 4: Generate separate MobiDBLite-style linear schematics

The matrix plots remain shared-domain dot matrices. To inspect where domains and
MobiDBLite disorder regions occur along each protein sequence, generate a
separate linear schematic and separate tables:

```bash
python3 scripts/plot_mobidblite_linear_schematic.py \
  --proteins data/output/polymerase_tf_proteins.csv \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --out data/output/all48_mobidblite_linear_schematic.svg \
  --regions-csv data/output/all48_mobidblite_linear_regions.csv \
  --summary-csv data/output/all48_mobidblite_linear_summary.csv \
  --summary-md data/output/all48_mobidblite_linear_summary.md
```

Outputs:
- `data/output/all48_mobidblite_linear_schematic.svg`
- `data/output/all48_mobidblite_linear_regions.csv`
- `data/output/all48_mobidblite_linear_summary.csv`
- `data/output/all48_mobidblite_linear_summary.md`

This schematic uses the existing protein sequences for protein lengths, raw
domain coordinates for domain blocks, and `MobiDBLite` source rows for disorder
blocks. Known SGD PTM sites are marked on the protein backbone when
`polymerase_tf_ptm_sites.csv` exists next to the `--proteins` file (or when
`--ptm-sites` points to one). Pass `--ptm-sites ""` to leave PTMs off. PTMs are
never added to the matrix plots.

### Step 5: Build the HTML output index

After the figures are written, build one HTML page that collects everything in
the output folder:

```bash
python3 scripts/build_output_index.py --outdir data/output
```

This writes `data/output/index.html`. Open it by double-clicking it (any
browser, macOS or Windows; no server or internet needed). It shows:

- a run summary from `run_summary.json`
- motif dot-matrix previews grouped into All proteins / Pol I / Pol II / Pol III
- buttons that open each figure's SVG, PNG, PDF, and order CSV
- each file's full path on your machine, with a "Copy path" button
- every table in the folder with a short description, row and column counts
- warnings for failed SGD/InterPro lookups or proteins missing from the tree

Links in the page are relative, so keep `index.html` inside the output folder.
The folder can then be moved, zipped, or shared between macOS and Windows.
Rerun the command whenever you add plots. Useful options:

- `--embed-previews`: embed the preview images in the HTML so the page still
  shows the figures if it is sent on its own (the file links still need the folder).
- `--recursive`: include files in subfolders (`cache/` is always skipped).
- `--html-out PATH` and `--title TEXT`: choose the file name and page title.

Each plotting command also writes a small `<figure>.plot.json` file next to the
SVG. It records the category, options (deduped, empty columns hidden, threshold),
inputs, and counts, and the index uses it to label each figure. Figures made
before this file existed are labeled from their file names.

## Optional MTPrompt-PTM Candidate Layer

MTPrompt-PTM is supported as an optional hypothesis-generation layer only. This
repository does not include MTPrompt source code, model weights, or datasets.

Prepare per-type FASTA inputs for the 13 supported MTPrompt types:

```bash
python3 scripts/run_mtprompt_ptm.py \
  --proteins data/output/polymerase_tf_proteins.csv \
  --outdir data/output/mtprompt_ptm \
  --prepare-only
```

Normalize externally generated prediction CSVs:

```bash
python3 scripts/run_mtprompt_ptm.py \
  --proteins data/output/polymerase_tf_proteins.csv \
  --known-ptms data/output/polymerase_tf_ptm_sites.csv \
  --prediction-dir data/output/mtprompt_ptm/predictions
```

Outputs:
- `data/output/polymerase_tf_ptm_predictions_mtprompt.csv`
- `data/output/polymerase_tf_ptm_candidate_sites.csv`

Candidate evidence tiers are kept separate:
- `known_sgd`: SGD known PTM site, whether or not MTPrompt predicts it.
- `predicted_known_overlap`: MTPrompt prediction at a matching SGD known site.
- `predicted_novel`: MTPrompt prediction with no matching SGD known site.

Predicted-only sites are never merged into known PTM counts.

## Category-Specific Plotting

You can render Pol I / Pol II / Pol III subsets from the same combined tree
(the script prunes the tree automatically). For exact category-specific
recreations, pass the matching category-specific tree if available.

### Example: Pol II only

```bash
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
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
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
  --style-config configs/plot_layout.compact.example.json \
  --out data/output/all48_motif_figure_compact.svg
```

Generate a fresh default template:

```bash
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
  --write-default-style-config configs/my_style.json
```

### Option B: Direct CLI overrides

```bash
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
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
- `--dedupe-harmonized`: collapse mapped motif/domain rows to InterPro IDs where possible; unmapped source motifs are retained as grey rows
- `--hide-empty-columns`: remove protein columns with no plotted motif/domain hits
- `--png`: optional high-resolution PNG export path. The Matplotlib renderer writes it directly; the native renderer rasterizes the SVG with `rsvg-convert` when available and otherwise falls back to Pillow
- `--pdf`: optional vector PDF export path. The native renderer uses `rsvg-convert` when available and otherwise needs ReportLab
- `--export-dpi`: PNG export DPI; default is 600
- Margins: `--margin-left/right/top/bottom`
- Cell geometry: `--cell-width`, `--cell-height`
- Vertical spacing: `--title-subtitle-gap`, `--subtitle-matrix-gap`, `--dendrogram-gap`
- Typography: `--font-family`, `--title-font-size`, `--subtitle-font-size`, `--protein-label-font-size`, `--domain-label-font-size`
- Colors: `--dot-fill`, `--dot-stroke`, `--unmapped-motif-color`, `--grid-color`, `--background-color`, `--category-color`
- Legend: `--no-legend`, `--legend-position`, `--legend-x`, `--legend-y`

For category-specific images, rerun the same plotting command with
`--category "Pol I"`, `--category "Pol II"`, or `--category "Pol III"` and set
category-specific `--out`, `--png`, and `--order-csv` paths. The renderer keeps
the plot structure intact when a category has no shared motifs at the selected
threshold, and the auto legend position is kept below the title/subtitle area
for small plots.

## Matplotlib Renderer

The Matplotlib renderer is the recommended plotting path. It uses the same
input tables, tree files, filtering, and most layout/style options as the
original native SVG renderer.

Additional notes are in `docs/matplotlib_renderer.md`.

```bash
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
  --tree data/input/clustalo-all48.phylotree \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --proteins data/output/polymerase_tf_proteins.csv \
  --out data/output/all48_motif_figure_matplotlib.svg \
  --png data/output/all48_motif_figure_matplotlib.png \
  --pdf data/output/all48_motif_figure_matplotlib.pdf
```

## Output and Caching Notes

- All generated tables/figures are intended to live in `data/output/`.
- `index.html` (Step 5) summarizes an output folder; `*.plot.json` files describe each figure.
- API responses are cached under `data/output/cache/`.
- `.gitignore` excludes generated outputs and common large assets from version control.

## Caveats

- SGD/InterPro responses can change over time, so exact row counts may drift.
- If your workbook schema differs from the expected columns, extraction will fail until the mapping is updated.
- Large generated images/PDFs should not be committed.
