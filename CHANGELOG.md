# Changelog

Notable changes to the pipeline, newest first. Each entry says what changed,
what it affects, and how to get the previous behavior when that is possible.

## 2026-09-28: one-command runs, HTML index, motif-based ordering, Pol I subgroups

### Added

- **One-command pipeline** (`scripts/run_pipeline.py`). It goes from a protein
  CSV to extraction tables, then dot matrices for all proteins and each Pol
  category, then the linear schematic, and finally `index.html`. The command
  is the same on macOS and Windows.
  - Each dot matrix is drawn raw and InterPro-deduped, in both column orders.
  - Files are named `<prefix>_<scope>_<raw|deduped>_by_<tree|motifs>`.
- **HTML output index** (`scripts/build_output_index.py`). Writes one
  `index.html` per output folder.
  - Figure previews are grouped into All proteins / Pol I / Pol II / Pol III.
  - Each figure links its SVG, PNG, PDF and column-order CSV, and shows its
    column order and suggested legend.
  - Full file paths are shown with copy buttons.
  - Every table is listed with a description and row/column counts.
  - Failed lookups and proteins missing from the tree are flagged.
  - It needs only the standard library and uses relative links, so it works on
    macOS and Windows and after the folder is moved or zipped.
- **Order columns by shared motifs** (`--order-by motifs` in both plotting
  scripts). This is average-linkage clustering on the Jaccard distance between
  the motif sets actually plotted.
  - Proteins with the same motifs sit together.
  - Proteins without plotted motifs go last.
  - The result does not depend on input row order.
  - No tree file is needed.
- **Figure manifests.** Both plotting scripts write `<figure>.plot.json` with
  category, options, column order, suggested legend, inputs, outputs and counts.
- **Local Clustal Omega matches EBI.** `run_clustalo_tree.py --backend local`
  now writes the same neighbour-joining "phylotree" and percent identity matrix
  that the EBI service returns. In testing, alignment and tree were
  byte-identical to EBI's (see `docs/clustalo_ebi_reproducibility.md`).
  `--local-tree guide` / `--ebi-tree-result guide` give Clustal Omega's guide
  tree instead. `--max-wait` stops polling a stuck EBI job (default 1 hour).
- **Pol I axial / peri-axial subgroups.** A `pol_i_subgroup` input column (with
  aliases) is carried through all extraction tables and counted in
  `run_summary.json`. Plots color labels, dots and legend entries as
  "Pol I axial" / "Pol I peri-axial". Order CSVs gain `pol_i_subgroup` and
  `display_group` columns.
- `scripts/build_pol1_extended_input.py`: merges a base protein CSV with a Pol I
  workbook, labels axial/peri-axial rows, and fills missing sequences from UniProt.
- `scripts/plot_mobidblite_linear_schematic.py`: a separate linear schematic of
  domains and MobiDBLite disorder along each protein, plus region and summary
  tables.
- New plotting options:
  - `--dedupe-harmonized`: one row per InterPro entry; unmapped source motifs
    are kept as grey rows.
  - `--hide-empty-columns`
  - `--unmapped-motif-color`
- `--category-color` now accepts any group name, e.g. `"Pol I peri-axial=#6baed6"`.
- README: quick start, Windows setup and command notes, column-order guidance.
  Also this changelog and an MIT `LICENSE`.

### Changed (affects how existing figures look)

- **Default title and subtitle.** The title is now "All proteins" (or
  "<category> proteins"). The subtitle now states the column order, e.g.
  "Ordered by Clustal Omega sequence similarity (neighbour-joining)"; it was
  "Hierarchical Clustering". Pass `--title` / `--subtitle` to override.
- **Dot colors follow each protein's group** (Pol I/II/III, axial/peri-axial)
  instead of a single green. For the old look, pass
  `--dot-fill "#178b1d" --dot-stroke "#0c4c10"`.
- The line under the title is no longer drawn by default
  (`title_line_width` is now 0).
- The default `--export-dpi` for PNGs is 600 (was 300).
- The left margin widens automatically for long motif labels. Very sparse
  category plots keep a minimum matrix and figure width.
- The auto legend stays in the top corner unless it would overlap the title or
  subtitle; then it moves below the subtitle.
- PTM overlays are no longer drawn on matrix plots. `--show-ptm-overlay` and
  `--ptm-intersections` are hidden and ignored. PTM tables are still written by
  the extractor.
- Native renderer: PNG and PDF are converted from the final SVG with
  `rsvg-convert` when it is available, so they match the SVG exactly.
- Grey (unmapped) motif rows are sorted by the ID shown on the axis.
- `run_clustalo_tree.py` never silently substitutes the guide tree when EBI's
  phylogenetic tree is missing.

### Fixed

- The linear schematic crashed while writing the regions CSV
  (`harmonized_key` column), leaving an empty regions table and no summaries.
- Domain hits past the end of the protein sequence are skipped. They used to be
  clamped into fake 1-aa regions.
- Linear schematic PTM options:
  - `--ptm-sites ""` now omits PTMs; it used to fail.
  - By default the PTM table is read from the same folder as `--proteins`.
  - Methylation sites get their own color instead of "other".
- Native renderer `--pdf` failed unless ReportLab was installed. It now writes
  a vector PDF with `rsvg-convert` first.
- Plotting and summary scripts read and write files as UTF-8. Windows defaults
  to cp1252, and Excel adds a byte-order mark to CSV files.
- `build_pol1_extended_input.py` used UniProt's first search hit without
  checking that it was the requested gene; the search also matches aliases.
  - A hit whose main gene name matches is now preferred.
  - An alias-only match is accepted only when unambiguous, and is recorded in
    `source_note`.
  - Otherwise the script stops.

### Known limitations

- Clustal Omega results depend on the order of sequences in the FASTA. The
  script writes them in CSV row order, so keep that order fixed. For proteins
  from different families the tree is weak; see
  `docs/clustalo_ebi_reproducibility.md`.
- Only tested on macOS so far; Windows handling (paths, encodings, the HTML
  index) is designed in but not yet run on a Windows machine.

## Earlier history

- 2026-05: PTM integration pipeline (SGD known PTMs, optional MTPrompt-PTM
  candidate layer); Matplotlib renderer made headless-safe and documented.
- 2026-04: Matplotlib renderer added and tuned to match the native SVG
  layout; optional PNG/PDF exports; CSV protein-list input; initial repository
  scaffold.
