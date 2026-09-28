# Clustal Omega: EBI Web Service vs Local Runs

The Clustal Omega web form is:

- https://www.ebi.ac.uk/jdispatcher/msa/clustalo

`scripts/run_clustalo_tree.py --backend ebi` submits the same FASTA through the
EMBL-EBI Job Dispatcher REST service that backs that form:

- https://www.ebi.ac.uk/Tools/services/rest/clustalo

Example command (EBI asks for a contact email address):

```bash
python3 scripts/run_clustalo_tree.py \
  --backend ebi \
  --email xxs410@case.edu \
  --proteins data/input/proteins.csv \
  --fasta data/input/clustalo.fasta \
  --alignment-out data/input/clustalo.ebi.clu \
  --out data/input/clustalo.ebi.phylotree \
  --distmat-out data/input/clustalo.ebi.distmat \
  --job-id-out data/input/clustalo.ebi-job.txt
```

Parameters sent in the default `--ebi-mode full-distance`:

- `stype=protein`
- `outfmt=clustal_num`
- `order=aligned`
- `guidetreeout=true`
- `dismatout=true`
- `mbed=false`
- `mbediteration=false`
- `iterations=0`
- `gtiterations=-1`
- `hmmiterations=-1`

The web form's default keeps mBed on. EBI only returns the identity matrix when
mBed is off, which is why the script turns it off by default. Use
`--ebi-mode web-default` to keep the web form's setting.

Downloaded EBI results:

- alignment: `aln-clustal_num`
- tree used for plotting: `phylotree` (the web page's "Phylogenetic Tree")
- percent identity matrix: `pim`

EBI also returns `tree` ("Guide Tree"), Clustal Omega's own guide tree. Its
leaves are all about the same distance from the root, so it looks nearly flat
under the dot matrix. Use `--ebi-tree-result guide` only if you want it.

## Local runs give the same result

`--backend local` (the default) runs Clustal Omega on your machine with EBI's
settings (`--full --output-order=tree-order --outfmt=clu --resno`). It then
computes the tree and identity matrix the way EBI does:

- EBI's `phylotree` is a neighbour-joining tree that ClustalW 2.1 computes
  from the finished alignment.
- The distance between two sequences is their share of mismatches over the
  columns where neither has a gap. There is no correction for multiple
  substitutions.
- Pairs that share no such column get distance 0.
- Pairs are joined using Saitou and Nei's selection criterion.

The following was verified on 2026-09-28 on a lab dataset of yeast transcription-related proteins:

- Local Clustal Omega 1.2.4 (the version EBI runs) produced an alignment
  byte-identical to EBI's.
- The local tree was byte-identical to EBI's `phylotree`.
- The local identity matrix had the same values as EBI's `pim`.
- mBed on and mBed off gave the same alignment for that set and for a
  second, smaller set.

`tests/test_run_clustalo_tree.py` repeats the tree comparison whenever the lab
example files and `clustalo` are present.

## Two cautions when using the tree

**Sequence order in the FASTA changes the result.** Clustal Omega is not
order-independent. The script writes the FASTA in the protein CSV's row order,
so a fixed CSV gives a fixed tree. The same sequences pasted into the web form
in a different order give a different tree. On the lab dataset,
reshuffled orders gave trees barely closer to the original than random trees.
The comparisons below ignore branch flips and root position; a flipped copy of
the tree scores 100% on each:

- About 3% of internal branches were shared.
- Four-protein agreement was 39.5%, against 33.7% for random trees. This test
  asks whether two trees group each set of four proteins the same way, and it
  is not thrown off by a single protein moving.
- Leaving out the 20 least stable proteins raised agreement only to 51%.

**Weakly related proteins give a weak tree.** In that dataset:

- 96% of protein pairs shared under 20% identity (median 11%).
- 7% of pairs shared no gap-free aligned column. The tree step treats those
  pairs as identical, which shows up as paired positive and negative branch
  lengths.
- Only one grouping, a pair of paralogs, was stable across reshuffled orders.

For such sets the dendrogram is a reproducible display order, not evidence of
evolutionary relationships. The plotting scripts label tree-ordered figures
accordingly. They can also order columns by shared motifs instead
(`--order-by motifs`), which does not depend on input order.
