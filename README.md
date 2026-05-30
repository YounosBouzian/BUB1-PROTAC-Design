# BUB1-PROTAC-Design

> **Computational design of the first selective BUB1 kinase PROTAC degraders**  
> Dr. Younos BOUZIAN · Medicinal Chemistry & AI Drug Discovery · 2026

---

## Overview

This repository documents a complete **computational PROTAC design campaign** targeting BUB1 kinase — a key regulator of the spindle assembly checkpoint (SAC) overexpressed in multiple human cancers.

Starting from **BAY-1816032** (IC₅₀ = 8 nM vs BUB1, Kd = 3.3 nM), the most potent and selective BUB1 inhibitor published, three PROTAC candidates were designed, built in RDKit, and profiled computationally before synthesis.

**No selective BUB1 PROTAC exists in PROTAC-DB** — this work represents a first-in-class design effort.

---

## Scientific Background

### Why Target BUB1?

| Property | Value |
|---|---|
| Target | BUB1 kinase (Budding Uninhibited by Benzimidazoles 1) |
| Function | Spindle Assembly Checkpoint (SAC) kinase; phosphorylates H2A-T120; recruits Sgo1/2 |
| Cancer relevance | Overexpressed in breast, lung, colorectal, and other solid tumours |
| Clinical indication | Triple-negative breast cancer (TNBC) — validated by BAY-1816032 + olaparib synergy data |
| PROTAC rationale | Catalytic inhibition alone insufficient — BUB1 scaffolding functions persist |

### Why PROTAC over Inhibitor?

```
BAY-1816032 (inhibitor):          BUB1-PROTAC:
  ✅ Blocks kinase activity          ✅ Destroys ALL BUB1 functions
  ❌ Leaves scaffold intact          ✅ Eliminates scaffolding
  ❌ Requires continuous dosing      ✅ Catalytic — substoichiometric dose
  ❌ Resistance via overexpression   ✅ Resistant cells degrade more protein
```

### Warhead — BAY-1816032

- **Scaffold:** Indazole-pyrimidine-methoxypyridine  
- **IC₅₀ vs BUB1:** 8 nM (biochemical); 29 nM (cellular H2A-pT120)  
- **Selectivity:** >100-fold over BUB1B and >390 other kinases  
- **Key feature:** Solvent-exposed hydroxyethoxy chain (-OCH₂CH₂OH) confirmed as exit vector by crystal structure PDB: **6F7B**  
- **CAS:** 1891087-61-8

---

## Repository Structure

```
BUB1-PROTAC-Design/
│
├── README.md                          ← You are here
│
├── notebooks/
│   └── BUB1_PROTAC_design.py         ← Main design + analysis script
│
├── data/
│   └── protac_candidates.csv         ← SMILES + properties of all candidates
│
├── results/
│   └── admet_summary.csv             ← SwissADME + pkCSM results
│
└── figures/
    ├── building_blocks.png            ← Warhead + E3 ligands
    ├── protac_candidates.png          ← All 3 PROTAC structures
    └── warhead_exit_vector.svg        ← BAY-1816032 exit vector highlighted
```

---

## PROTAC Candidates

Three CRBN-recruiting PROTACs designed from BAY-1816032 via N-alkylation of the exit vector hydroxyl:

| ID | Linker | E3 Ligase | MW (Da) | LogP | TPSA (Å²) | Synthesis Priority |
|---|---|---|---|---|---|---|
| **PROTAC-1** | PEG2 + amide | CRBN (Pomalidomide) | 861.8 | 4.80 | 218.1 | 3rd |
| **PROTAC-2** | PEG3 + amide | CRBN (Pomalidomide) | 905.9 | 4.82 | 227.3 | **1st** |
| **PROTAC-3** | PEG4 + amide | CRBN (Pomalidomide) | 949.9 | 4.83 | 236.5 | 2nd |

### SMILES

```
# BAY-1816032 (warhead)
COc1cnccc1Nc1nc(-c2nn(Cc3c(F)cc(OCCO)cc3F)c3ccccc23)ncc1OC

# PROTAC-1 (PEG2/CRBN)
COc1cnccc1Nc1nc(-c2nn(Cc3c(F)cc(OCCOCCC(=O)Nc4cccc5c4C(=O)N(C4CCC(=O)NC4=O)C5=O)cc3F)c3ccccc23)ncc1OC

# PROTAC-2 (PEG3/CRBN) — Priority candidate
COc1cnccc1Nc1nc(-c2nn(Cc3c(F)cc(OCCOCCOCCC(=O)Nc4cccc5c4C(=O)N(C4CCC(=O)NC4=O)C5=O)cc3F)c3ccccc23)ncc1OC

# PROTAC-3 (PEG4/CRBN)
COc1cnccc1Nc1nc(-c2nn(Cc3c(F)cc(OCCOCCOCCOCCC(=O)Nc4cccc5c4C(=O)N(C4CCC(=O)NC4=O)C5=O)cc3F)c3ccccc23)ncc1OC
```

---

## Computational Results

### Physicochemical Properties (RDKit)

| Property | PROTAC-1 | PROTAC-2 | PROTAC-3 | PROTAC Rule |
|---|---|---|---|---|
| MW (Da) | 861.8 ✅ | 905.9 ✅ | 949.9 ✅ | < 1100 |
| LogP | 4.80 ✅ | 4.82 ✅ | 4.83 ✅ | 1–6 |
| TPSA (Å²) | 218.1 ✅ | 227.3 ✅ | 236.5 ✅ | < 250 |
| HBD | 3 ✅ | 3 ✅ | 3 ✅ | ≤ 5 |
| HBA | 14 ✅ | 15 ✅ | 16 ✅ | 5–20 |
| Rot. bonds | 16 ✅ | 19 ✅ | 22 ✅ | 10–25 |
| Fsp3 | 0.23 ✅ | 0.27 ✅ | 0.30 ✅ | > 0.20 |
| Tanimoto vs warhead | 0.546 ✅ | 0.546 ✅ | 0.546 ✅ | 0.25–0.75 |

### ADMET Profile (SwissADME + pkCSM)

| Property | PROTAC-1 | PROTAC-2 | PROTAC-3 |
|---|---|---|---|
| Aqueous solubility | 1.97e-4 mg/mL | 2.25e-4 mg/mL | 2.56e-4 mg/mL |
| GI absorption* | Low | Low | Low |
| Pgp substrate | Yes ⚠️ | Yes ⚠️ | Yes ⚠️ |
| CYP3A4 inhibitor | Yes ⚠️ | Yes ⚠️ | **No ✅** |
| Caco-2 (log cm/s) | 1.131 ✅ | 1.105 ✅ | 1.079 ✅ |
| hERG inhibition | No ✅ | No ✅ | No ✅ |
| Hepatotoxicity | No ✅ | No ✅ | No ✅ |
| AMES toxicity | No ✅ | No ✅ | No ✅ |
| Rule-based score | 1.00 ✅ | 1.00 ✅ | 1.00 ✅ |

*GI absorption prediction unreliable for PROTACs — tools trained on MW < 500 Da compounds. See ARV-110 (34% oral bioavailability in clinic despite "Low" prediction).

---

## Synthesis Plan

All three candidates are accessible in **2 steps from commercial materials**:

```
Step 1 — O-alkylation
  BAY-1816032 (-OH)  +  Br-(CH₂CH₂O)n-COOH
  K₂CO₃, DMF, 60°C, 12h
  → BAY-1816032-PEGn-COOH intermediate

Step 2 — Amide coupling
  BAY-1816032-PEGn-COOH  +  H₂N-Pomalidomide
  EDC·HCl / HOBt, DIPEA, DMF, RT, 12h
  → PROTAC product

Commercial sources:
  BAY-1816032:      MedKoo Biosciences / Sigma-Aldrich
  Pomalidomide:     Sigma-Aldrich / TCI
  PEG linkers:      Quanta Biodesign / BroadPharm
```

---

## How to Run

```bash
# Clone repository
git clone https://github.com/YounosBouzian/BUB1-PROTAC-Design.git
cd BUB1-PROTAC-Design

# Install dependencies
pip install rdkit pandas

# Run design pipeline
python notebooks/BUB1_PROTAC_design.py
```

**Or run directly in Google Colab:**  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com)

```python
!pip install rdkit --quiet
# Then paste BUB1_PROTAC_design.py content
```

---

## Key References

1. **Siemeister et al.** *Clin. Cancer Res.* 2019 — BAY-1816032 characterisation and BUB1 biology
2. **Bemis et al.** *J. Med. Chem.* 2021 — Physicochemical property analysis of 1,200 PROTACs
3. **Crew et al.** *J. Med. Chem.* 2018 — BRD4 PROTAC design methodology
4. **Bondeson et al.** *Nat. Chem. Biol.* 2015 — Catalytic PROTAC degradation concept
5. **Vamathevan et al.** *Nat. Rev. Drug Discov.* 2019 — AI in drug discovery overview

---

## Author

**Dr. Younos BOUZIAN**  
Medicinal Chemistry · AI-Assisted Drug Discovery · Python for Cheminformatics  

- 8+ years synthetic & SAR chemistry · 260+ compounds synthesised  
- Drug discovery pipeline: hit identification → Phase 1 clinical trials  
- Patent: PCT/EP2023/055381  

[![GitHub](https://img.shields.io/badge/GitHub-YounosBouzian-black?logo=github)](https://github.com/YounosBouzian)

---

## License

MIT License — academic and research use. If you use this work, please cite this repository.

---

*Generated as part of a 6-month AI + Drug Discovery self-development programme · May 2026*
