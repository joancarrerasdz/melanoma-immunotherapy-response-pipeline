# Biological interpretation of the strict 12-gene signature

## Purpose and scope

This document provides the biological interpretation of the 12-gene
consensus core produced by the strict leakage-controlled modelling
pipeline:

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

The purpose is interpretation, not post-hoc model optimization.

No gene was added, removed, reweighted, or selected on the basis of the
literature review reported here. The signature was defined by the
statistical pipeline before this biological interpretation.

The evidence is intentionally separated into:

1. established molecular or cellular function;
2. relevant melanoma or immune-checkpoint context;
3. hypotheses that remain speculative.

This signature must be considered a **candidate research signature**.
It is not a clinically validated biomarker and no causal relationship
with anti-PD-1 response is claimed.

## Main biological signal: B-cell-associated biology

Four genes in the strict signature form the clearest coherent biological
group:

- **FCRL1**
- **FCRL2**
- **MS4A1**
- **PLCG2**

FCRL1 and FCRL2 are Fc-receptor-like proteins expressed predominantly in
B-cell compartments and participate in regulation of B-cell receptor
(BCR) signalling. Experimental work supports a signal-enhancing role for
FCRL1, whereas FCRL2 contains inhibitory motifs and can suppress
BCR-mediated signalling.

MS4A1 encodes CD20, a canonical B-cell surface molecule involved in
B-cell development and differentiation.

PLCG2 encodes phospholipase C gamma 2, an intracellular signalling enzyme
that generates IP3 and DAG downstream of immune-system receptors,
including BCR-associated signalling pathways.

Taken together, the presence of these four genes in a stable model-derived
gene set is biologically compatible with a contribution from B-cell
abundance, activation state, differentiation state, or spatial immune
organization in the tumour microenvironment.

This interpretation is consistent with independent melanoma literature.
Studies by Helmink et al. and Cabrita et al. reported enrichment of
B-cell-associated transcriptional programs and tertiary lymphoid
structures (TLS) in melanomas with favourable immune-checkpoint-blockade
outcomes.

However, the present data are bulk transcriptomic data. They cannot
determine whether the four genes reflect:

- greater numbers of infiltrating B cells;
- a specific B-cell functional state;
- tertiary lymphoid structures;
- another correlated immune-cell program; or
- a combination of these factors.

Therefore, the strict signature does **not** demonstrate that B cells or
TLS mechanistically cause anti-PD-1 response in these cohorts.

## Innate immune context: TLR1

**TLR1** encodes Toll-like receptor 1, a pattern-recognition receptor with
an established role in innate immune recognition and cytokine-inducing
signalling.

Its inclusion is compatible with variation in innate immune activity or
immune-cell composition in the tumour microenvironment.

The current analysis does not establish which cell population contributes
the TLR1 signal and does not demonstrate a direct role for TLR1 in
anti-PD-1 sensitivity.

Accordingly, TLR1 is interpreted as an immune-context candidate rather
than a validated response mechanism.

## Proliferation and DNA replication context: TICRR

**TICRR** (Treslin) participates in the initiation of DNA replication
through its interaction with replication-control machinery.

Experimental cancer studies have associated increased TICRR activity with
enhanced DNA replication and tumour-cell proliferation in multiple solid
tumour contexts.

Within this signature, TICRR may therefore reflect variation in
proliferative or replication-related cellular states.

This interpretation is not specific to melanoma and no evidence from the
current analysis demonstrates that TICRR directly modifies anti-PD-1
response.

## Post-transcriptional and cell-state context: TRIM71

**TRIM71** encodes an E3 ubiquitin ligase and RNA-binding regulatory
protein involved in post-transcriptional regulation and cell-cycle-related
processes.

Its presence may represent variation in cellular state or proliferative
programs.

There is insufficient evidence in the present analysis to connect TRIM71
directly to anti-PD-1 response in melanoma. Any such connection should
therefore be considered hypothesis-generating only.

## Genes with limited response-specific interpretability

The remaining genes have identifiable molecular functions but currently
lack sufficiently strong evidence to assign them a specific melanoma
anti-PD-1 mechanism.

### DENND2A

DENND2A encodes a DENN-domain guanine-nucleotide exchange factor involved
in intracellular trafficking, including endosome-to-Golgi transport.

This provides a valid cellular-function annotation, but a direct
connection to anti-PD-1 response is not established here.

### RASL10B

RASL10B belongs to the Ras-like small GTPase family.

Its response-specific functional interpretation remains limited, and the
current study does not support assigning a particular immune or melanoma
mechanism to this gene.

### SCARF2

SCARF2 encodes a scavenger-receptor-family protein with reported roles
related to receptor interactions and cell adhesion.

Evidence is insufficient to interpret SCARF2 as a specific determinant of
anti-PD-1 response.

### SDR42E1

SDR42E1 is a short-chain dehydrogenase/reductase family member with
predicted oxidoreductase activity and predicted involvement in steroid
biosynthetic processes.

These annotations do not provide sufficient evidence for a specific
melanoma immunotherapy mechanism.

### VARS2

VARS2 encodes the mitochondrial valyl-tRNA synthetase required for
mitochondrial translation.

Its inclusion could reflect variation in mitochondrial or general cellular
state, but the current study provides no evidence that this constitutes a
specific anti-PD-1 response mechanism.

## Integrated interpretation

The strongest interpretable pattern in the strict 12-gene set is the
co-occurrence of multiple genes linked to B-cell identity and BCR
regulation:

**FCRL1, FCRL2, MS4A1 and PLCG2.**

This pattern is biologically compatible with published observations that
B-cell-rich immune microenvironments and tertiary lymphoid structures can
be associated with immune-checkpoint-blockade response in melanoma.

A second, weaker layer comprises genes connected to broader immune or
cell-state biology, including TLR1, TICRR and TRIM71.

The remaining genes should be retained as statistically selected members
of the strict signature without forcing a mechanistic interpretation that
is not supported by available evidence.

Importantly, the model was derived from bulk RNA sequencing. Consequently,
gene-expression signals may encode tumour-cell biology, immune-cell
abundance, stromal composition, cellular state, or correlated mixtures of
these processes.

## What can be claimed

The results support the following statement:

> The strict model-derived signature contains a coherent B-cell-associated
> component, alongside genes linked to innate immune signalling,
> replication/cell state and other cellular processes. This biological
> structure is compatible with independent evidence implicating B-cell-rich
> tumour microenvironments in melanoma checkpoint-blockade response.

## What cannot be claimed

The analysis does not establish that:

- any of the 12 genes is individually predictive in independent clinical
  practice;
- any gene is causally responsible for anti-PD-1 response;
- the signature directly measures tertiary lymphoid structures;
- the signature is independent of immune-cell composition;
- the 12-gene panel is a validated clinical biomarker;
- targeting any of these genes would improve immunotherapy response.

External performance remains materially lower than internal performance,
so biological plausibility must not be interpreted as clinical validation.

## Evidence classification

| Gene | Primary context | Evidence strength for functional annotation | Anti-PD-1 interpretation |
|---|---|---:|---|
| FCRL1 | B-cell/BCR regulation | Strong | Plausible immune-context link; not causal |
| FCRL2 | B-cell/BCR regulation | Strong | Plausible immune-context link; not causal |
| MS4A1 | B-cell identity/CD20 | Strong | Strong contextual relevance to B-cell-rich melanoma microenvironments |
| PLCG2 | Immune-receptor signalling | Strong | Plausible BCR/immune-context link |
| TLR1 | Innate immune recognition | Strong | Plausible immune-context link |
| TICRR | DNA replication/proliferation | Strong | General tumour/cell-state context; not response-specific |
| TRIM71 | Post-transcriptional/cell-state regulation | Moderate | Hypothesis-generating only |
| DENND2A | Vesicular trafficking | Moderate | Response-specific evidence insufficient |
| RASL10B | Small-GTPase biology | Limited | Response-specific evidence insufficient |
| SCARF2 | Scavenger receptor/cell adhesion | Moderate | Response-specific evidence insufficient |
| SDR42E1 | Oxidoreductase/metabolic biology | Limited/predicted | Response-specific evidence insufficient |
| VARS2 | Mitochondrial translation | Strong molecular function | Response-specific evidence insufficient |

## References

1. Helmink BA, et al. B cells and tertiary lymphoid structures promote
   immunotherapy response. Nature. 2020;577:549-555.
   doi:10.1038/s41586-019-1922-8. PMID: 31942075.

2. Cabrita R, et al. Tertiary lymphoid structures improve immunotherapy
   and survival in melanoma. Nature. 2020;577:561-565.
   doi:10.1038/s41586-019-1914-8. PMID: 31942071.

3. Zhao X, et al. Fc receptor-like 1 intrinsically recruits c-Abl to
   enhance B cell activation and function. Science Advances. 2019.
   doi:10.1126/sciadv.aaw0315. PMID: 31328160.

4. Jackson TA, et al. FcR-like 2 inhibition of B cell receptor-mediated
   activation of B cells. Journal of Immunology. 2010.
   doi:10.4049/jimmunol.1002305. PMID: 21068405.

5. NCBI Gene 931. MS4A1 — membrane spanning 4-domains A1.

6. NCBI Gene 5336. PLCG2 — phospholipase C gamma 2.

7. NCBI Gene 7096. TLR1 — toll like receptor 1.

8. NCBI Gene 90381. TICRR — TOPBP1 interacting checkpoint and
   replication regulator.

9. Yu Q, et al. TICRR contributes to tumorigenesis through accelerating
   DNA replication in cancers. Frontiers in Oncology. 2019.
   PMID: 31275851.

10. NCBI Gene 27147. DENND2A — DENN domain containing 2A.

11. NCBI Gene 91179. SCARF2 — scavenger receptor class F member 2.

12. NCBI Gene 93517. SDR42E1 — short chain dehydrogenase/reductase
    family 42E member 1.

13. NCBI Gene 131405. TRIM71 — tripartite motif containing 71.

14. NCBI Gene 57176. VARS2 — valyl-tRNA synthetase 2, mitochondrial.

## QA interpretation

This document addresses portfolio QA criterion **P1-11**:

> Genes connected with biological mechanisms through literature,
> pathways or functional evidence, with explicit separation between
> evidence and speculation.

The interpretation is downstream of the frozen statistical analysis and
does not modify any model, feature-selection decision, threshold, metric,
or released analytical result.
