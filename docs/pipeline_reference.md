# Pipeline Reference

## Scripts

- `scripts/extract_sgd_domains.py`
  - Reads workbook + calls SGD/InterPro APIs
  - Produces protein/domain tables, matrices, source summary, and cache files

- `scripts/summarize_shared_domains.py`
  - Aggregates shared raw domains by category composition

- `scripts/plot_clustalo_motif_figure.py`
  - Builds SVG motif matrix + dendrogram from tree/hits/proteins tables
  - Supports category-specific plots and layout customization via CLI/JSON

## Typical Data Flow

1. Workbook (`data/input/*.xlsx`) -> extraction outputs in `data/output/`
2. Raw hits table -> shared-domain summary table
3. Tree + raw hits + proteins table -> SVG figure + order CSV
