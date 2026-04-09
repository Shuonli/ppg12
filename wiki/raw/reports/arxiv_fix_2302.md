# arXiv ID Fix: 2302.01560 in Proposed Reading List

## Problem

The file `wiki/proposed_reading_list.md` (line 62, section A2 LHC) contains the entry:

> | 2 | 2302.01560 | Isolated photon pp 13 TeV, 36 fb-1 (JHEP 07, 2023) | ATLAS | pp 13 TeV |

This arXiv ID is **wrong**, and the journal reference is also incorrect.

## What arXiv:2302.01560 actually is

- **Title:** "Describe, Explain, Plan and Select: Interactive Planning with Large Language Models Enables Open-World Multi-Task Agents"
- **Authors:** Zihao Wang, Shaofei Cai, Guanzhou Chen, Anji Liu, Xiaojian Ma, Yitao Liang
- **Subject:** Artificial Intelligence (cs.AI)
- **Conference:** NeurIPS 2023
- **DOI:** 10.48550/arXiv.2302.01560

This is a computer science paper about LLM-based planning agents in Minecraft. It has no connection to ATLAS, photon physics, or particle physics.

## Correct arXiv ID

The intended paper is **arXiv:1908.02746**.

- **Title:** Measurement of the inclusive isolated-photon cross section in pp collisions at sqrt(s) = 13 TeV using 36 fb^{-1} of ATLAS data
- **Collaboration:** ATLAS
- **Journal reference:** JHEP 10 (2019) 203
- **DOI:** 10.1007/JHEP10(2019)203
- **arXiv category:** hep-ex
- **Submitted:** 2019-08-07
- **INSPIRE record:** 1772071
- **INSPIRE URL:** https://inspirehep.net/literature/1772071

## What needs to be corrected in the reading list

| Field | Wrong (current) | Correct |
|-------|----------------|---------|
| arXiv ID | 2302.01560 | 1908.02746 |
| Journal ref | JHEP 07 (2023) | JHEP 10 (2019) 203 |

The corrected row should read:

> | 2 | 1908.02746 | Isolated photon pp 13 TeV, 36 fb-1 (JHEP 10, 2019) | ATLAS | pp 13 TeV | Latest precision ATLAS photon; state-of-art NLO comparison technique |

## Verification method

1. Fetched https://arxiv.org/abs/2302.01560 -- confirmed it is an AI/NeurIPS paper
2. Queried INSPIRE-HEP API: `find t isolated photon and cn atlas and j JHEP and de 2022-2024` (with fields filter for titles, arxiv_eprints, publication_info, dois)
3. Identified arXiv:1908.02746 as the unique match for "ATLAS inclusive isolated-photon cross section, 13 TeV, 36 fb-1, JHEP"
4. Cross-verified via arXiv API (`export.arxiv.org/api/query?id_list=1908.02746`) -- title, journal ref, and DOI all confirmed

## Note on "JHEP 07 2023" claim

No ATLAS isolated-photon paper was published in JHEP 07 (2023). The reading list entry appears to have both a fabricated arXiv ID and an incorrect journal reference. The actual publication date is October 2019, not July 2023.
