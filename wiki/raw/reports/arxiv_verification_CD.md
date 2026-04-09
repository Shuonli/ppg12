# arXiv Verification Report: Sections C and D

Verified 2026-04-08 using `arxiv.org/abs/{ID}` (title, authors, journal ref) and `arxiv.org/html/{ID}` (HTML availability).

## Tier Definitions

- **Tier 1**: HTML available at `arxiv.org/html/{ID}` -- full-text ingestible
- **Tier 2**: Abstract page only at `arxiv.org/abs/{ID}` -- metadata + abstract ingestible
- **Tier 3**: No arXiv ID -- manual entry required

---

## Section C: Theory & Methods -- Isolation, NLO, Fragmentation, PDFs

| arXiv ID | Verified Title | Authors/Collab | Journal Ref | Correct? | Tier | Notes |
|----------|----------------|----------------|-------------|----------|------|-------|
| hep-ph/0602133 | A new critical study of photon production in hadronic collisions | P. Aurenche, M. Fontannaz, J.-Ph. Guillet, E. Pilon, M. Werlen | Phys.Rev. D73 (2006) 094007 | Yes | 2 | JETPHOX NLO program. Reading list says "PRD 73, 2006" -- matches. |
| hep-ph/9704447 | Quarks and Gluon Fragmentation Functions into Photons | L. Bourhis, M. Fontannaz, J.Ph. Guillet | Eur.Phys.J. C2 (1998) 529-537 | Yes | 2 | BFG fragmentation functions. Reading list says "EPJC 2, 1998" -- matches. |
| hep-ph/9801442 | Isolated photons in perturbative QCD | S. Frixione | Phys.Lett. B429 (1998) 369-374 | Yes | 2 | Smooth-cone isolation. Reading list says "PLB 429, 1998" -- matches. |
| hep-ph/0204023 | Cross section of isolated prompt photons in hadron-hadron collisions | S. Catani, M. Fontannaz, J.Ph. Guillet, E. Pilon | JHEP 0205 (2002) 028 | Yes | 2 | NLO isolated photon framework. Reading list says "JHEP 05, 2002" -- matches. |
| 1506.07443 | New parton distribution functions from a global analysis of quantum chromodynamics | S. Dulat, T.J. Hou, J. Gao, M. Guzzi, J. Huston, P. Nadolsky, J. Pumplin, C. Schmidt, D. Stump, C.P. Yuan | Phys.Rev. D93 (2016) 033006 | Yes | 2 | CT14 PDFs. Reading list says "PRD 93, 2016" -- matches. |
| 1005.4529 | Sensitivity of isolated photon production at TeV hadron colliders to the gluon distribution in the proton | R. Ichou, D. d'Enterria | Phys.Rev. D82 (2010) 014015 | Yes | 2 | Cone vs smooth isolation. Reading list says "PRD 82, 2010" -- matches. |
| 2312.02519 | **Creative Agents: Empowering Agents with Imagination for Creative Tasks** | P. Cai, C. Zhang, Y. Fu, H. Yuan, Z. Lu | PMLR 244:471-496 (UAI 2025) | **WRONG ID** | N/A | **ERROR**: This is an AI/ML paper, not the NNLO isolated photon paper. See notes below. |
| Gordon:1993qc | (no arXiv) | L.E. Gordon, W. Vogelsang | PRD 48 (1993) | N/A | 3 | Vogelsang NLO calculation. No arXiv ID (pre-arXiv era). |
| RMP 59 (1987) 465 | (no arXiv) | J.F. Owens | Rev.Mod.Phys. 59 (1987) 465 | N/A | 3 | Classic prompt photon review. No arXiv ID (pre-arXiv era). |

### Issue: 2312.02519 is incorrect

The reading list claims this is "NNLO isolated photon (PRL 132, 2024)" by "Campbell, Ellis, Seth". Verification reveals:

1. **arXiv 2312.02519** is actually "Creative Agents: Empowering Agents with Imagination for Creative Tasks" by Cai et al. -- a computer science paper completely unrelated to QCD.

2. **The closest matching paper found** is arXiv **1612.04333**: "Direct Photon Production at Next-to-Next-to-Leading Order" by J.M. Campbell, R.K. Ellis, C. Williams, published in **PRL 118 (2017) 222001**. This is the first NNLO calculation for direct photon production.

3. **No paper matching "PRL 132, 2024" by Campbell, Ellis, Seth** was found in INSPIRE-HEP. The coauthor may be "Williams" not "Seth", and the journal ref may be "PRL 118, 2017" not "PRL 132, 2024".

4. There are also NNLO isolated photon papers by Chen, Gehrmann, Glover, Hofer, Huss: arXiv **1904.01044** (JHEP 04, 2020) and **1905.08577** (2019 proceedings). These are the "isolated" (not just "direct") photon NNLO calculations.

**Recommended action**: The user should verify on INSPIRE which specific NNLO photon paper they intend to cite. Candidates:
- **1612.04333** -- Campbell, Ellis, Williams, PRL 118 (2017) 222001 [direct photon at NNLO]
- **1904.01044** -- Chen, Gehrmann, Glover, Hofer, Huss, JHEP 04 (2020) 166 [isolated photon at NNLO]

---

## Section D: Techniques -- Unfolding, MVA/BDT, Background Subtraction

| arXiv ID | Verified Title | Authors/Collab | Journal Ref | Correct? | Tier | Notes |
|----------|----------------|----------------|-------------|----------|------|-------|
| 1010.0632 | Improved iterative Bayesian unfolding | G. D'Agostini | (unpublished, conference proceedings) | Yes | 2 | Reading list says "D'Agostini improved Bayesian unfolding (2010)" -- matches. No formal journal ref; presented at Alliance Workshop on Unfolding, Hamburg 2010. |
| 1105.1160 | Unfolding algorithms and tests using RooUnfold | T. Adye | CERN-2011-006, pp 313-318 (PHYSTAT 2011) | Yes | 2 | Reading list says "RooUnfold (PHYSTAT 2011)" -- matches. Conference proceedings. |
| physics/0703039 | TMVA -- Toolkit for Multivariate Data Analysis | A. Hoecker et al. (27 authors) | CERN-OPEN-2007-007 | Yes | 2 | Reading list says "TMVA toolkit" -- matches. Technical report, not a journal paper. |
| 1603.02754 | XGBoost: A Scalable Tree Boosting System | T. Chen, C. Guestrin | KDD 2016 (DOI: 10.1145/2939672.2939785) | Yes | 2 | Reading list says "XGBoost (KDD 2016)" -- matches. |
| 1012.4389 | Measurement of the inclusive isolated prompt photon cross section in pp collisions at sqrt(s) = 7 TeV with the ATLAS detector | ATLAS Collaboration | Phys.Rev. D83 (2011) 052005 | Yes | 2 | Reading list says "ATLAS first 7 TeV photon with ABCD (PRD 83, 2011)" -- matches. Pioneered 2x2 ABCD sideband method. |
| hep-ph/9509307 | SVD Approach to Data Unfolding | A. Hoecker, V. Kartvelishvili | Nucl.Instrum.Meth. A372 (1996) 469-481 | Yes | 2 | Reading list says "SVD unfolding" -- matches. |
| NIM A 362 (1995) 487 | (no arXiv) | G. D'Agostini | Nucl.Instrum.Meth. A362 (1995) 487 | N/A | 3 | Original Bayesian unfolding paper. Pre-arXiv. |

---

## Summary

| Section | Total Papers | Verified OK | Wrong ID | No arXiv (Tier 3) | All Tier 2 (abs only) |
|---------|-------------|-------------|----------|--------------------|-----------------------|
| C | 9 | 6 | 1 (2312.02519) | 2 | Yes (all with arXiv IDs) |
| D | 7 | 6 | 0 | 1 | Yes (all with arXiv IDs) |
| **Total** | **16** | **12** | **1** | **3** | -- |

### HTML Availability (Tier 1 vs Tier 2)

All 12 papers with arXiv IDs were checked at `arxiv.org/html/{ID}` -- every one returned HTTP 404. This is expected: arXiv HTML conversion (ar5iv/HTML5) primarily covers papers submitted after ~2023, and all Section C/D papers predate this. **All arXiv papers are Tier 2 (abstract + PDF only).**

### Action Items

1. **CRITICAL**: Replace arXiv ID **2312.02519** with the correct NNLO photon paper. Most likely candidate is **1612.04333** (Campbell, Ellis, Williams, PRL 118, 2017) or **1904.01044** (Chen et al., JHEP 04, 2020). Verify which paper is intended.
2. All other 12 arXiv IDs verified correct against their reading list descriptions.
3. Three papers have no arXiv ID (Tier 3): Gordon:1993qc, RMP 59 (1987) 465, NIM A 362 (1995) 487. These require manual BibTeX entries.
