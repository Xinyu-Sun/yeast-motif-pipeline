# Matplotlib Renderer

This is the recommended plotting implementation for the
motif/domain matrix figure:

- Script: `scripts/plot_clustalo_motif_figure_matplotlib.py`
- Library: Matplotlib
- Inputs: tree, hits, proteins, and optional style config files
- Outputs: SVG by default, with optional high-resolution PNG/vector PDF via
  `--png` and `--pdf`
- PTM extraction tables are kept separate and are not rendered on matrix plots
- If the proteins table contains `pol_i_subgroup`, Pol I labels and legend
  entries distinguish axial and peri-axial groups.
- `--dedupe-harmonized` collapses mapped rows to InterPro motif/domain IDs
  where possible and keeps unmapped source motifs as grey rows.
- `--hide-empty-columns` removes protein columns with no plotted motif/domain
  hits after the selected filters are applied.
- Long y-axis motif labels expand the left margin automatically so deduped
  InterPro labels are not clipped.
- The auto legend stays in the top corner unless it would overlap the centered
  title or subtitle (small or sparse plots); then it moves below the subtitle.
- Each run also writes `<figure>.plot.json`, which `scripts/build_output_index.py`
  uses to label the figure in the HTML output index.

Install the plotting dependency:

```bash
python3 -m pip install -r requirements-plotting.txt
```

Example:

```bash
python3 scripts/plot_clustalo_motif_figure_matplotlib.py \
  --tree data/input/clustalo-all48.phylotree \
  --hits data/output/polymerase_tf_domain_hits_raw.csv \
  --proteins data/output/polymerase_tf_proteins.csv \
  --out data/output/all48_motif_figure_matplotlib.svg \
  --png data/output/all48_motif_figure_matplotlib.png \
  --pdf data/output/all48_motif_figure_matplotlib.pdf \
  --order-csv data/output/all48_motif_order_matplotlib.csv \
  --min-proteins 2
```

Notes:

- It preserves the same ordering and domain filtering behavior as the native
  SVG renderer.
- For exact Pol I / Pol II / Pol III recreations, use the matching
  category-specific tree when available.
