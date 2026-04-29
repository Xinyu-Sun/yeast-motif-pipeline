# Pipeline Reference

## Scripts

- `scripts/extract_sgd_domains.py`
  - Reads a CSV protein table or legacy workbook + calls SGD/InterPro APIs
  - Produces protein/domain tables, matrices, source summary, and cache files

- `scripts/summarize_shared_domains.py`
  - Aggregates shared raw domains by category composition

- `scripts/plot_clustalo_motif_figure.py`
  - Builds SVG motif matrix + dendrogram from tree/hits/proteins tables
  - Optionally exports PNG and PDF copies of the plot
  - Supports category-specific plots and layout customization via CLI/JSON

- `scripts/plot_clustalo_motif_figure_matplotlib.py`
  - Experimental Matplotlib implementation of the same figure
  - Reuses the same inputs, ordering, filtering, and style controls

## Typical Data Flow

1. Protein list (`data/input/*.csv`) or workbook (`data/input/*.xlsx`) -> extraction outputs in `data/output/`
2. Raw hits table -> shared-domain summary table
3. Tree + raw hits + proteins table -> SVG figure, optional PNG/PDF exports, and order CSV
