# Reproducible environment contract

## Purpose

This document formalizes the reproducible R environment for the strict
Week 3 nested-cross-validation and Week 4 locked external-validation
workflow.

The environment snapshot was created during QA remediation after Week 4
had already been completed. At snapshot time, the installed versions of
the core modelling packages exactly matched the environment recorded
when the frozen Week 4 deployment model was created.

No Week 3 or Week 4 analytical result was recomputed or modified in order
to create this environment contract.

## Historical Week 4 frozen environment

The original Week 4 deployment freeze records:

- R 4.5.1
- edgeR 4.6.3
- caret 7.0.1
- randomForest 4.7.1.2
- pROC 1.19.0.1

Source:

`results/week4_deployment_environment_strict.csv`

These versions are treated as historical frozen provenance and are
checked against `renv.lock`.

## Direct dependencies of the strict W3-W4 workflow

The explicit project dependency contract is declared in `DESCRIPTION`.

Direct dependencies are:

- AnnotationDbi
- caret
- edgeR
- org.Hs.eg.db
- pROC
- randomForest
- readxl

Transitive dependencies are captured recursively in `renv.lock`.

The project uses an explicit renv snapshot so that unrelated historical
or exploratory scripts do not expand the dependency contract.

## Bioconductor

The project is locked to Bioconductor 3.21.

This is recorded by the renv project configuration and lockfile.

## Restore

From the repository root, start with the R version recorded in
`renv.lock` and run:

    Rscript -e 'renv::restore(prompt = FALSE)'

The project `.Rprofile` activates renv for ordinary R sessions.

After restoration, the strict workflows can be executed with:

    Rscript scripts/run_week3.R
    Rscript scripts/run_week4.R

The historical `Rscript --vanilla ...` commands remain valid records of
how the workflows were previously executed, but `--vanilla` bypasses
project startup files and therefore should not be used when the purpose
is specifically to execute inside the restored renv environment.

## Verification

Run:

    Rscript scripts/19_reproducible_environment_gate.R

The gate verifies:

1. required renv infrastructure exists;
2. the dependency contract uses explicit snapshot mode;
3. all required direct dependencies are declared;
4. the lockfile records Bioconductor 3.21;
5. the historical Week 4 R/package versions match the lockfile;
6. the project-local renv package library is not tracked by Git.

A successful check terminates with:

    U11 REPRODUCIBLE ENVIRONMENT: PASS

## Methodological boundary

The reproducibility snapshot does not redefine the analytical workflow.

In particular, it does not:

- change preprocessing;
- reselect genes;
- rerun model development to obtain different results;
- retune Week 4 hyperparameters;
- alter the frozen classification threshold;
- replace the historical Week 4 environment record.

Its role is to make the existing strict W3-W4 implementation
reconstructible and auditable.
