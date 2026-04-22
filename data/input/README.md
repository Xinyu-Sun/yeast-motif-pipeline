# Input Files

Place raw user-provided inputs in this folder.

Required for the full pipeline:
- `polymerase_tf_sequences.xlsx`: source workbook containing Pol I/II/III + combined protein blocks.
- `clustalo-all48.phylotree`: Clustal Omega guide tree in Newick format.

Included sample:
- `clustalo-all48.phylotree` (copied from the original figure recreation workspace).

Notes:
- Workbook files are ignored by `.gitignore` to avoid committing local or sensitive data.
- The plotting script can also use a different tree path via `--tree`.
