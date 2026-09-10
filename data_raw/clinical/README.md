# GSE160638 clinical data source

This directory records the authoritative patient-level clinical source used to reconstruct the GSE160638 endpoint.

## Official supplementary source

- Title: *Table S2 from MYC Induces Immunotherapy and IFNγ Resistance Through Downregulation of JAK2*
- Repository: AACR Figshare
- Publisher: American Association for Cancer Research
- Persistent DOI: <https://doi.org/10.1158/2326-6066.30721063>
- Landing page: <https://aacr.figshare.com/articles/dataset/Table_S2_from_MYC_Induces_Immunotherapy_and_IFN_Resistance_Through_Downregulation_of_JAK2/30721063>
- Accessed on: 2026-09-01

## Associated publication

- Markovits E, Harush O, Baruch EN, et al.
- Title: *MYC Induces Immunotherapy and IFNγ Resistance Through Downregulation of JAK2*
- Journal: *Cancer Immunology Research*
- Year: 2023
- Volume and pages: 11(7):909–924
- DOI: <https://doi.org/10.1158/2326-6066.CIR-22-0184>

## Role in this project

Table S2 is the authoritative source for patient-level clinical outcomes in GSE160638. The response endpoint must not be inferred from free-text GEO metadata.

Only worksheet `S2A` will be used for the primary anti-PD-1 cohort. The download URL, local filename, file size, SHA-256 checksum, and source structure will be recorded after completing the automated-download step.

## Reuse note

The workbook is publicly accessible through AACR Figshare. No specific reuse licence is inferred in this repository; users must consult the source record and associated publication for applicable terms.

## Week 1 analytical use

For the Week 1 rebuild, worksheet `S2A` is parsed programmatically by `scripts/01_load_data.R`.

The clinical source is used to:

- define the primary anti–PD-1 population;
- preserve the official `CR/PR/SD/PD` response endpoint;
- derive the binary `Responder/NonResponder` endpoint;
- normalize and validate PD1 patient identifiers;
- align clinical records with GEO metadata and expression counts.

The validated primary cohort contains 36 unique PD1 samples:

- 22 Responders;
- 14 NonResponders.

Samples belonging to the TIL-ACT/TIL component of the source dataset (`TIL_*`) are not included in the primary anti–PD-1 cohort.

Clinical missingness and category inconsistencies are preserved and documented rather than silently corrected during Week 1.

For the complete analytical-data specification, see the repository-level `DATA.md`.
