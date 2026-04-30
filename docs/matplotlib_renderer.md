# Matplotlib Renderer

This renderer provides a standard plotting-library implementation for the
motif/domain matrix figure:

- Script: `scripts/plot_clustalo_motif_figure_matplotlib.py`
- Library: Matplotlib
- Inputs: same tree, hits, proteins, and style config files as the native SVG
  renderer
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

- This is intentionally additive: it does not replace the native SVG renderer.
- It is useful when you want Matplotlib-managed SVG, PNG, or PDF output while
  preserving the same ordering and domain filtering as the native renderer.
