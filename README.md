# clinvar-sas-reclassification
A bioinformatics pipeline for ancestry-aware interpretation of ClinVar hereditary cancer variants using gnomAD South Asian population data.

# 🧬 Reclassification of ClinVar Hereditary Cancer Variants Using South Asian Population Data

## Overview

This repository contains the complete computational workflow developed as part of an M.Sc. Bioinformatics dissertation.

The project investigates how South Asian population-specific allele frequency data from gnomAD v4.1 can improve hereditary cancer variant interpretation using ACMG/AMP population-frequency criteria and ClinGen Variant Curation Expert Panel (VCEP) recommendations.

---

## Dissertation Title

**Reclassification of ClinVar Hereditary Cancer Variants Using South Asian Population Frequency Data from gnomAD v4.1**

---

## Research Objectives

- Evaluate the utility of South Asian (SAS) population frequency data.
- Integrate ClinVar and gnomAD v4.1 datasets.
- Apply BA1, BS1 and PM2_Supporting ACMG/AMP criteria.
- Resolve Variants of Uncertain Significance (VUS).
- Resolve conflicting ClinVar interpretations.
- Identify pathogenic-frequency discordant variants.
- Develop a transparent and reproducible bioinformatics workflow.

---

## Genes Analysed

- APC
- ATM
- BRCA1
- BRCA2
- CDH1
- MLH1
- MSH2
- PALB2
- PTEN
- TP53

---

## Workflow

ClinVar Records
↓
Variant Parsing & Normalization
↓
MySQL Database Construction
↓
gnomAD v4.1 SAS Matching
↓
FAF95 Calculation
↓
BA1 / BS1 / PM2_Supporting Assignment
↓
Variant Reclassification
↓
Statistical Analysis
↓
Visualization & Reporting

---

## Technologies Used

- Python
- MySQL
- Pandas
- NumPy
- SciPy
- Matplotlib

---

## Key Results

| Metric | Value |
|----------|----------|
| Total Variants Analysed | 102,825 |
| Found in gnomAD SAS | 32,620 |
| Absent from gnomAD SAS | 70,205 |
| BA1 Assigned | 279 |
| BS1 Assigned | 820 |
| PM2_Supporting Assigned | 69,068 |
| Genuine Reclassifications | 572 |
| VUS Resolved | 126 |
| Conflicting Variants Resolved | 345 |
| False Positive Candidates | 5 |

---

## Repository Structure

project/
├── scripts/
├── results/
├── figures/
├── thesis/
├── README.md
├── requirements.txt

---

## Installation

```bash
git clone <repository-url>
cd repository
pip install -r requirements.txt
```

---


---

## Author

Prince Kushwaha

M.Sc. Bioinformatics

Guru Gobind Singh Indraprastha University (GGSIPU)

New Delhi, India
