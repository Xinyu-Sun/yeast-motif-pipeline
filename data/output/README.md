# Output Files

Pipeline-generated outputs are written here by default.

Expected outputs:
- `polymerase_tf_proteins.csv`
- `polymerase_tf_domain_hits_raw.csv`
- `polymerase_tf_domain_hits.csv`
- `polymerase_tf_domain_hits_harmonized.csv`
- `polymerase_tf_domain_source_summary.csv`
- `polymerase_tf_domain_matrix_*.csv`
- `all48_shared_domain_category_summary.csv`
- figure outputs such as `*.svg` and order tables such as `*_order.csv`
- `*.plot.json` figure descriptions written by the plotting scripts
- `index.html` from `scripts/build_output_index.py`: open it in a browser to browse all figures and tables
- cached API responses in `cache/`

This directory is intentionally git-ignored (except this README) to keep the repository lightweight.
