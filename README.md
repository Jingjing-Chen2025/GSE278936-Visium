# GSE278936-Visium

Spatial Transcriptomic analysis on GSE278936 Visium — Full 7-Tier Multi-Sample Pipeline

## Overview

This repository contains a complete R pipeline for analyzing 52 Visium spatial transcriptomics samples from GSE278936, covering:

- **TIER 1:** Spatial Annotation (Define the Interface)
- **TIER 2:** Distance-Binned Pseudo-Bulk DE (Cross-Patient Statistics)
- **TIER 3:** Gradient Modeling (Gene-by-Gene Decay Analysis)
- **TIER 4:** Ligand-Receptor Communication (Paracrine Inference)
- **TIER 5:** Trajectory / Transition State (Spatial Pseudotime)
- **TIER 6:** Gene Regulatory Network (Interface-Specific GRNs)
- **TIER 7:** Cross-Sample Meta-Analysis (Conserved Interface Programs)

## Signatures

- ASPC (Adipose Stromal Progenitor Cell)
- Hallmark AR (Androgen Response)
- NEPC Beltran UP (Neuroendocrine Prostate Cancer upregulated)
- NEPC Beltran DOWN (Neuroendocrine Prostate Cancer downregulated)

## Samples

52 samples across 5 conditions: BPH (4), TRNA (17), NEADT (22), CRPC (5), MET (4)