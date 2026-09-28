# Pipeline Reference

## Scripts

- `scripts/extract_sgd_domains.py`
  - Reads a CSV protein table or legacy workbook + calls SGD/InterPro APIs
  - Produces protein/domain tables, matrices, source summary, and cache files

- `scripts/build_pol1_extended_input.py`
  - Merges the existing 48-protein CSV with the Pol I peri-axial workbook
  - Labels original Pol I rows as axial and newly added workbook rows as peri-axial
  - Fills missing sequences from UniProt with cached lookups

- `scripts/summarize_shared_domains.py`
  - Aggregates shared raw domains by category composition

- `scripts/run_clustalo_tree.py`
  - Writes FASTA from a protein CSV and runs local or EMBL-EBI Clustal Omega
  - Produces the alignment, EBI's neighbour-joining "phylotree", and a percent identity matrix; local and EBI outputs are identical for the same FASTA
  - See `docs/clustalo_ebi_reproducibility.md` for the EBI parameters and the local-vs-EBI comparison

- `scripts/plot_clustalo_motif_figure.py`
  - Native SVG implementation of the motif matrix + dendrogram
  - Orders columns by the Clustal Omega tree (`--order-by tree`) or by average-linkage clustering on shared motifs (`--order-by motifs`)
  - Optionally exports PNG and PDF copies of the plot
  - Supports category-specific plots and layout customization via CLI/JSON
  - Supports category-colored dots, InterPro-deduped motif rows, and hiding empty protein columns

- `scripts/plot_clustalo_motif_figure_matplotlib.py`
  - Recommended Matplotlib implementation of the same figure
  - Reuses the same inputs, ordering, filtering, and style controls
  - Supports the same InterPro-deduped and no-empty-column matrix options

- `scripts/run_pipeline.py`
  - Runs extraction, shared-domain summary, all-protein and per-category dot matrices (raw and InterPro-deduped), the linear schematic, and the HTML index in one command

- `scripts/build_output_index.py`
  - Writes a standalone `index.html` for an output folder: figure previews grouped by All / Pol I / Pol II / Pol III, links to SVG/PNG/PDF/CSV files, full paths with copy buttons, table descriptions, and warnings
  - Uses relative links so the folder works after moving or sharing between macOS and Windows

- `scripts/plot_mobidblite_linear_schematic.py`
  - Produces a separate protein-length schematic for domain coordinates and MobiDBLite disorder regions
  - Writes independent SVG, region CSV, per-protein summary CSV, and Markdown summary outputs
  - Does not alter motif/domain matrix plots and does not add PTM information to them

## Typical Data Flow

1. Optional Pol I workbook merge -> extended protein CSV with subgroup labels and sequences
2. Protein list (`data/input/*.csv`) or workbook (`data/input/*.xlsx`) -> extraction outputs in `data/output/`
3. Protein table -> Clustal Omega FASTA, guide tree, alignment, and distance matrix
4. Raw hits table -> shared-domain summary table
5. Tree + raw hits + proteins table -> Matplotlib SVG figure, optional PNG/PDF exports, and order CSV
6. Proteins table + raw hits table -> separate MobiDBLite linear schematic outputs and summaries
7. Output folder -> `index.html` summarizing all figures and tables
