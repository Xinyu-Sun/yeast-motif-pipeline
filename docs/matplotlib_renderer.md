# Matplotlib Renderer

This is the recommended plotting implementation for the
motif/domain matrix figure:

- Script: `scripts/plot_clustalo_motif_figure_matplotlib.py`
- Library: Matplotlib
- Inputs: tree, hits, proteins, and optional style config files
- Outputs: SVG by default, with optional PNG/PDF via `--png` and `--pdf`

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
