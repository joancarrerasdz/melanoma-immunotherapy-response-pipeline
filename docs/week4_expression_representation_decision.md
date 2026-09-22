# Week 4 — External-validation expression representation decision

## Objective

Week 4 requires a common predictive representation for the frozen strict
Week 3 core genes across:

- GSE160638 training cohort;
- GSE91061 external cohort;
- GSE78220 external cohort.

The training cohort originates from RNA-seq raw counts and uses TMM/logCPM
internally, whereas the external cohorts provide FPKM expression values.

Direct numerical equivalence between these expression scales must therefore
not be assumed.

## Frozen gene space

The external-validation gene space is restricted to the 12 genes selected in
all five strict Week 3 outer folds:

- DENND2A
- FCRL1
- FCRL2
- MS4A1
- PLCG2
- RASL10B
- SCARF2
- SDR42E1
- TICRR
- TLR1
- TRIM71
- VARS2

All 12 genes are available in both external cohorts.

## Candidate representations

Three label-free representations were audited:

1. native log scale
   - GSE160638: TMM/logCPM;
   - external cohorts: log2(FPKM + 1);

2. sample-wise rank transformation;

3. sample-wise z-score transformation.

No external response labels were used during this comparison.

## Cross-dataset similarity

### GSE91061 versus GSE160638

| Representation | Pearson | Spearman | Mean absolute gene-median difference |
|---|---:|---:|---:|
| native_log_scale | 0.7718 | 0.8252 | 1.5329 |
| samplewise_rank | 0.8589 | 0.8891 | 0.0990 |
| samplewise_z | 0.7854 | 0.8671 | 0.3226 |

### GSE78220 versus GSE160638

| Representation | Pearson | Spearman | Mean absolute gene-median difference |
|---|---:|---:|---:|
| native_log_scale | 0.7678 | 0.8671 | 1.6196 |
| samplewise_rank | 0.9125 | 0.9330 | 0.0851 |
| samplewise_z | 0.7647 | 0.9021 | 0.3648 |

## Decision

`samplewise_rank` is selected as the Week 4 deployment representation.

The transformation is calculated independently within each sample using only
the frozen 12-gene strict Week 3 core.

This representation was selected because it:

- avoids assuming numerical equivalence between TMM/logCPM and FPKM;
- is invariant to monotonic differences in expression scale;
- requires no external response labels;
- preserves within-sample relative gene-expression information;
- showed the strongest training-to-external gene-profile concordance in both
  external cohorts.

## Methodological boundary

The Week 3 validation metrics remain associated exclusively with the strict
nested-cross-validation pipeline implemented in Week 3.

In particular, the Week 3 AUC of 0.8896 must not be presented as performance
of the 12-gene sample-wise-rank deployment model.

The 12-gene core was derived from the Week 3 training cohort and is treated as
a candidate deployment signature.

Any internal modelling performed after fixing this signature and
representation is part of model development and model freezing, not
independent validation.

Independent performance assessment will be performed only after the external
prediction pipeline, hyperparameters, probability threshold and prediction
outputs have been frozen without consulting external response labels.

## Gate result

**WEEK 4 EXPRESSION REPRESENTATION DECISION: PASS**

Frozen representation: `samplewise_rank`

External response labels used for representation selection: **FALSE**

External performance calculated during representation selection: **FALSE**
