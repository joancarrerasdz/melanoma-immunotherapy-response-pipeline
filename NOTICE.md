# Licensing and third-party material

## Original software

Original source code authored for this repository is released under the
MIT License. See `LICENSE`.

This primarily applies to the project-authored analytical and QA code in
`scripts/`.

## Third-party datasets

The MIT License does not grant rights over datasets obtained from external
sources.

In particular, the repository uses or derives analytical inputs from public
research resources including:

- GSE160638;
- GSE91061;
- GSE78220;
- associated supplementary clinical material.

These resources remain subject to the terms, attribution requirements and
conditions of their original providers and publications.

See `DATA.md` for provenance.

## Data and derived analytical artifacts

The following repository areas are not being relicensed as third-party data
under the MIT License:

- `data_raw/`
- `data_processed/`
- `results/`

Some files in these directories are derived from externally sourced research
data and are retained for provenance, reproducibility and portfolio review.

## Figures

Figures in `figures/` are project-generated analytical outputs. Their
underlying numerical information may derive from externally sourced datasets.

Their presence in this repository should not be interpreted as relicensing
the underlying source datasets.

## Historical thesis material

Historical master's-thesis material and preserved thesis outputs are separate
from the software licensing boundary.

This includes, where present:

- `docs/Memoria.pdf`;
- `docs/historical_thesis_results.md`;
- material preserved through the `tfm-original` tag.

## R packages and environment metadata

Third-party R packages retain their own respective licenses.

Files such as `renv.lock` record dependency and environment metadata and do
not relicense the packages they reference.
