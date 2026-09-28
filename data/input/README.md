# Input Files

Place raw user-provided inputs in this folder.

Required for the full pipeline:
- `polymerase_tf_sequences.xlsx`: source workbook containing Pol I/II/III + combined protein blocks.
- `clustalo-all48.phylotree`: Clustal Omega guide tree in Newick format.

Included sample:
- `clustalo-all48.phylotree` (copied from the original figure recreation workspace).
- `proteins.example.csv`: template for CSV-based protein-list input.

## CSV Protein Table

For a simpler input format, create a CSV table and pass it to the extractor with
`--csv`.

Required lookup column:
- Prefer `y_name`: SGD locus/systematic name used for API lookup.
- A single `protein`, `name`, or `standard_name` column also works if its values are valid SGD lookup names.

Optional columns:
- `category`: group label such as `Pol I`, `Pol II`, or `Pol III`.
- `protein`: display name.
- `accession`: accession identifier.
- `pol_i_subgroup`: optional Pol I subgroup such as `axial` or `periaxial`.
- `fasta_header`: FASTA header without sequence.
- `amino_acid_sequence`: protein sequence without FASTA header.
- `raw_fasta`: complete FASTA cell, including header and sequence.

Accepted aliases include `Y Name`, `systematic_name`, or `locus` for `y_name`;
`name` or `standard_name` for `protein`; and `sequence` or `aa_sequence` for
`amino_acid_sequence`.

Notes:
- Workbook files are ignored by `.gitignore` to avoid committing local or sensitive data.
- Local CSV input files are ignored by `.gitignore`; example files named
  `*.example.csv` are tracked.
- Tree files (`*.phylotree`, `*.nwk`) and EBI job-id files are ignored too,
  except the original `clustalo-all48.phylotree`.
- `scripts/build_pol1_extended_input.py` can generate an all-48 plus Pol I
  axial/peri-axial CSV from the Pol I peri-axial workbook.
- `scripts/run_clustalo_tree.py` can regenerate `.fasta`, `.clu`, `.distmat`,
  and `.phylotree` files from any protein CSV with sequences, using either a
  local Clustal Omega binary or the EMBL-EBI Job Dispatcher REST service.
- The plotting script can also use a different tree path via `--tree`.
