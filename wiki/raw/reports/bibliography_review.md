# PPG12 Bibliography Review

Generated: 2026-04-08

Reviewed files:
- `PPG12-analysis-note/cite.bib` (11 entries)
- `ppg12_conf_note/references.bib` (72 entries)
- `wiki/proposed_reading_list.md` (curated reading list, ~70 papers)

---

## 1. Must-Cite Gaps (three-star papers from the reading list missing from `cite.bib`)

The reading list identifies 12 three-star must-add papers. Here is the status of the 8 specifically flagged:

| Paper | arXiv / Ref | In `cite.bib`? | In `references.bib`? | Status |
|-------|-------------|----------------|----------------------|--------|
| TMVA toolkit | physics/0703039 | **NO** | **NO** | **MISSING from both** |
| XGBoost | 1603.02754 | **NO** | **NO** | **MISSING from both** |
| D'Agostini 1995 original | NIM A 362 (1995) 487 | **NO** | **NO** | **MISSING from both** |
| CT14 PDFs | 1506.07443 | **NO** | **NO** | **MISSING from both** |
| Frixione smooth cone isolation | hep-ph/9801442 | **NO** | **NO** | **MISSING from both** |
| Ichou & d'Enterria | 1005.4529 | **NO** | **NO** | **MISSING from both** |
| ATLAS ABCD (first 7 TeV photon) | 1012.4389 | **NO** | **NO** | **MISSING from both** |
| sPHENIX first paper | 2504.02240 | **NO** | **NO** | **MISSING from both** |

**All 8 flagged papers are confirmed missing from both bib files.** These are critical omissions:
- TMVA and XGBoost are directly used tools (TMVA::Experimental::RBDT in `apply_BDT.C`, XGBoost in `main_training.py`).
- CT14lo is explicitly used in the JETPHOX NLO predictions (mentioned in `results.tex` line 21).
- Frixione isolation is the theoretical foundation for the smooth-cone isolation used in NLO calculations compared to the data.
- D'Agostini 1995 is the original Bayesian unfolding paper; only the 2010 update is cited.
- The ATLAS ABCD paper pioneered the background subtraction method used in PPG12.
- The sPHENIX first paper serves as the de facto NIM for detector performance.

### Additional three-star papers also missing

| Paper | arXiv / Ref | Status |
|-------|-------------|--------|
| PHENIX first direct photon pp 200 GeV | hep-ex/0502006 | Missing from both |
| PHENIX isolated direct photon pp 200 GeV (only prior RHIC iso photon) | 1006.1347 | Missing from both |
| PHENIX direct photon ALL pp 510 GeV | 2202.08158 | Missing from both |
| NLO isolated photon (Catani et al.) | hep-ph/0204023 | Missing from both |

---

## 2. Conf Note to Analysis Note Sync

`references.bib` contains 72 entries. The conf note (`ppg12_conf_note/main.tex`) cites 31 unique keys. Of those, **25 keys are NOT in `cite.bib`**. These should be ported to the analysis note bib as the analysis note matures:

### Critical (used in theory comparison / world data plots)

| Key | Paper | Why port |
|-----|-------|----------|
| `Aurenche:2006vj` | JETPHOX NLO (hep-ph/0602133) | Theory comparison in results; `results.tex` references JETPHOX but only cites BFG, not JETPHOX itself |
| `Gordon:1993qc` | Vogelsang NLO (PRD 48, 1993) | Werner Vogelsang calculation mentioned in `results.tex` line 21; no citation for it |
| `ATLAS:2011ezy` | ATLAS isolated photon 7 TeV (35 pb-1) | World data comparison |
| `ATLAS:2013sdn` | ATLAS isolated photon 7 TeV (4.6 fb-1) | World data comparison |
| `ATLAS:2016fta` | ATLAS isolated photon 8 TeV | World data comparison |
| `ATLAS:2017nah` | ATLAS isolated photon 13 TeV | World data comparison |
| `CMS:2018qao` | CMS isolated photon + photon+jets 13 TeV | World data comparison |
| `CMS:2020oen` | CMS isolated photon PbPb+pp 5.02 TeV | World data comparison |
| `ALICE:2019rtd` | ALICE isolated photon 7 TeV | World data comparison |
| `ALICE:2024kgy` | ALICE isolated photon 13 TeV | World data comparison |
| `ALICE:2024yvg` | ALICE isolated photon pp+PbPb 5.02 TeV | World data comparison |
| `ALICE:2025bnc` | ALICE isolated photon pp+pPb | World data comparison |
| `ATLAS:2019ery` | ATLAS prompt photon p+Pb | World data comparison |
| `ATLAS:2015rlt` | ATLAS isolated photon PbPb 2.76 TeV | World data comparison |

### Detector and methods

| Key | Paper | Why port |
|-----|-------|----------|
| `Belmont:2023fau` | sPHENIX predictions (NPA 1043) | Physics case |
| `PHENIX:2015siv` | sPHENIX upgrade proposal (1501.06197) | Foundational document |
| `Sphenix:TDR` | sPHENIX TDR | Detector design |
| `HARRISON2003235` | RHIC project overview | Accelerator |
| `Aidala:2020toz` | EMCal beam test (IEEE TNS 68) | Detector validation |
| `sPHENIX:2017lqb` | EMCal+HCal beam test (IEEE TNS 65) | Detector validation |
| `Klest:2020sdb` | sPHENIX TPC | Detector subsystem |
| `Aune:2024ysr` | sPHENIX TPOT | Detector subsystem |
| `OConnor:1998llb` | BaBar solenoid magnet | Magnet reference |
| `David:2000wfa` | PHENIX PbSc pattern recognition | Calorimeter technique |
| `DAgostini:2010hil` | D'Agostini improved unfolding (alternate key) | Already covered by `dagostini2010improvediterativebayesianunfolding` |

**Highest priority for immediate porting:** `Aurenche:2006vj` (JETPHOX) and `Gordon:1993qc` (Vogelsang). The analysis note text in `results.tex` explicitly discusses these NLO calculations but provides no `\cite{}` for either.

---

## 3. Dangling Citations

Citations in `.tex` files with no corresponding entry in `cite.bib`:

| Key | File | Line | Status |
|-----|------|------|--------|
| `MBDxsec` | `systematics.tex` | 126 | **No bib entry in either file.** Will render as `[?]`. Needs the actual sPHENIX internal note / measurement for the MBD minimum-bias cross section. |

This is the only dangling citation. All other `\cite{}` keys in the analysis note (`Adye:2011gm`, `AGOSTINELLI2003250`, `Aguilar:2021sfa`, `Bourhis:1997yu`, `dagostini2010improvediterativebayesianunfolding`, `PHENIX:2012jgx`, `ppg09IAN`, `Sjostrand:2014zea`, `Spousta_2016`) resolve to entries in `cite.bib`.

---

## 4. Duplicate and Stale Entries

### Duplicate papers with different keys in `references.bib`

These are the same underlying paper appearing under two different bib keys:

| Paper | Key 1 | Key 2 | Resolution |
|-------|-------|-------|------------|
| D'Agostini improved unfolding (1010.0632) | `dagostini2010improvediterativebayesianunfolding` | `DAgostini:2010hil` | Keep one; conf note cites `DAgostini:2010hil`, analysis note cites the other |
| sPHENIX upgrade proposal (1501.06197) | `sPHENIX_proposal` | `PHENIX:2015siv` | Keep `PHENIX:2015siv` (used in conf note); remove `sPHENIX_proposal` |
| GEANT4 (NIM A 506) | `AGOSTINELLI2003250` | `GEANT4` | Keep `AGOSTINELLI2003250` (used in both notes); remove `GEANT4` |
| sPHENIX predictions (NPA 1043) | `sphenix_predictions` | `Belmont:2023fau` | Keep `Belmont:2023fau` (used in conf note); remove `sphenix_predictions` |
| EMCal beam test (IEEE TNS 68) | `EM_test_beam` | `Aidala:2020toz` | Keep `Aidala:2020toz` (used in conf note); remove `EM_test_beam` |
| EMCal+HCal beam test (1704.01461) | `EMCal_HCal_test_beam` | `sPHENIX:2017lqb` | Keep `sPHENIX:2017lqb` (used in conf note); `EMCal_HCal_test_beam` unused |

### Stale / uncited entries in `cite.bib`

| Key | Paper | Note |
|-----|-------|------|
| `Lis1` | sPHENIX Centrality Calibration | Not cited in any `.tex` file |
| `Seidlitz1` | sPHENIX EMCal Calibrations | Not cited in any `.tex` file |

These may be intentionally kept for future sections. If the analysis note is near final, they can be removed.

### Stale / uncited entries in `references.bib`

The following `references.bib` entries are not cited in the conf note `main.tex`:

`ADLER2001488`, `alice_detdeta`, `ampt`, `atlas_detdeta`, `bjorken_energy_density`, `Brahms_pi_spectra`, `Brahms_p_spectra`, `Busza_2018`, `cms_pbpb_detdeta`, `cms_pPb_detdeta`, `Cunqueiro_2022`, `dagostini2010improvediterativebayesianunfolding`, `EM_test_beam`, `EMCal_HCal_test_beam`, `epos`, `Eremin_2003`, `GEANT4`, `hijing`, `Karsch_2002`, `Kharzeev_2001`, `Lis1`, `Loizides_2015`, `PHENIX_200gev`, `PHENIX_detdeta`, `PHENIX_npq_scaling_detdeta`, `phenix_centrality`, `phenix_centrality_2`, `phenix_ebj_meas`, `phenix_mbd_nim`, `phenix_particle_spectra`, `Romatschke_2017`, `Seidlitz1`, `single_hadrons`, `sphenix_predictions`, `sPHENIX_proposal`, `STAR_detdeta`, `star_lambda_spectra`, `Tiwari_Sahoo_2017`, `Trainor_2015`, `Wang_Gyulassy_2001`

These 40 entries appear to be inherited from a different sPHENIX analysis (dET/deta) and are not relevant to the PPG12 isolated photon measurement. They can be pruned from `references.bib`.

---

## 5. Generated BibTeX Entries for Top 5 Most Critical Missing References

### 5.1. TMVA Toolkit (physics/0703039)

```bibtex
@article{Hocker:2007ht,
    author = "Hoecker, Andreas and Speckmayer, Peter and Stelzer, Joerg and Therhaag, Jan and von Toerne, Eckhard and Voss, Helge",
    title = "{TMVA: Toolkit for Multivariate Data Analysis}",
    eprint = "physics/0703039",
    archivePrefix = "arXiv",
    primaryClass = "physics.data-an",
    reportNumber = "CERN-OPEN-2007-007",
    journal = "PoS",
    volume = "ACAT",
    pages = "040",
    year = "2007"
}
```

### 5.2. XGBoost (1603.02754)

```bibtex
@inproceedings{Chen:2016btl,
    author = "Chen, Tianqi and Guestrin, Carlos",
    title = "{XGBoost: A Scalable Tree Boosting System}",
    booktitle = "{Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining}",
    eprint = "1603.02754",
    archivePrefix = "arXiv",
    primaryClass = "cs.LG",
    doi = "10.1145/2939672.2939785",
    pages = "785--794",
    year = "2016"
}
```

### 5.3. CT14 Parton Distributions (1506.07443)

```bibtex
@article{Dulat:2015mca,
    author = "Dulat, Sayipjamal and Hou, Tie-Jiun and Gao, Jun and Guzzi, Marco and Huston, Joey and Nadolsky, Pavel and Pumplin, Jon and Schmidt, Carl and Stump, Daniel and Yuan, C.-P.",
    title = "{New parton distribution functions from a global analysis of quantum chromodynamics}",
    eprint = "1506.07443",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "SMU-HEP-15-12, CTEQ-1507, MSUHEP-150707, PITT-PACC-1507",
    doi = "10.1103/PhysRevD.93.033006",
    journal = "Phys. Rev. D",
    volume = "93",
    number = "3",
    pages = "033006",
    year = "2016"
}
```

### 5.4. Frixione Smooth Cone Isolation (hep-ph/9801442)

```bibtex
@article{Frixione:1998jh,
    author = "Frixione, Stefano",
    title = "{Isolated photons in perturbative QCD}",
    eprint = "hep-ph/9801442",
    archivePrefix = "arXiv",
    reportNumber = "ETH-TH-98-01",
    doi = "10.1016/S0370-2693(98)00454-7",
    journal = "Phys. Lett. B",
    volume = "429",
    pages = "369--374",
    year = "1998"
}
```

### 5.5. D'Agostini 1995 Original Bayesian Unfolding

```bibtex
@article{DAgostini:1994fjx,
    author = "D'Agostini, G.",
    title = "{A multidimensional unfolding method based on Bayes' theorem}",
    doi = "10.1016/0168-9002(95)00274-X",
    journal = "Nucl. Instrum. Meth. A",
    volume = "362",
    pages = "487--498",
    year = "1995"
}
```

---

## Summary of Recommended Actions

### Immediate (before next draft)

1. **Add all 8 flagged three-star papers** to `cite.bib` (BibTeX entries for top 5 provided above).
2. **Fix `MBDxsec` dangling citation** in `systematics.tex` -- add the actual sPHENIX internal note for MBD cross-section to `cite.bib`.
3. **Port `Aurenche:2006vj` (JETPHOX) and `Gordon:1993qc` (Vogelsang)** from `references.bib` to `cite.bib` and add `\cite{}` calls in `results.tex` where the calculations are discussed.
4. **Add `\cite{}` for CT14lo** in `results.tex` line 21 where "CT14lo parton distribution functions" is mentioned.

### Short-term

5. **Port the 14 world-data measurement entries** from `references.bib` to `cite.bib` as comparison plot citations are added to the analysis note.
6. **Port detector references** (`Sphenix:TDR`, `PHENIX:2015siv`, `Belmont:2023fau`, `HARRISON2003235`, `Aidala:2020toz`, `sPHENIX:2017lqb`) as the detector section is written.
7. **Resolve duplicate bib keys** in `references.bib` (6 pairs identified above).

### Cleanup

8. **Prune ~40 stale entries** from `references.bib` inherited from dET/deta analysis.
9. **Remove or cite `Lis1` and `Seidlitz1`** in `cite.bib`.
