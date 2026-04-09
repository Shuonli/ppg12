---
title: "Pattern recognition in the PHENIX PbSc electromagnetic calorimeter"
authors: "G. David, E. Kistenev, S. White, C. Woody, A. Bazilevsky, S. Belikov, V. Kochetkov, V. Onuchin, A. Usachev"
journal: "IEEE Trans. Nucl. Sci. 47 (2000) 1982-1986"
category: E_detector
priority: 1
status: stub_needs_enrichment
tags: [PHENIX, EMCal, clustering, pattern-recognition, calorimeter, PbSc]
ppg12_relevance: "Documents EMCal clustering algorithms for tower-based calorimeters at RHIC, providing methodological heritage for the sPHENIX EMCal cluster reconstruction used in PPG12"
---

# Pattern Recognition in the PHENIX PbSc Electromagnetic Calorimeter

## Summary

This paper describes the pattern recognition and clustering algorithms developed for the PHENIX lead-scintillator (PbSc) electromagnetic calorimeter. The algorithms address the challenge of identifying and separating electromagnetic showers from individual photons and electrons in a segmented tower calorimeter, particularly in high-multiplicity environments. The methods include seed-tower identification, nearest-neighbor clustering, and shower-shape-based splitting of merged clusters. These algorithms established the framework for EMCal cluster reconstruction at RHIC that has influenced subsequent detector generations including sPHENIX.

## Key Points

- Clustering algorithm for segmented tower EMCal: seed tower above threshold, iterative addition of neighboring towers above a lower threshold
- Cluster splitting: when a cluster contains multiple local maxima, the shower is split by assigning tower energies to sub-clusters using expected shower profiles
- Shower shape characterization: cluster moments (width in eta and phi directions) used to distinguish single photons from merged pi0 decays and hadronic showers
- Performance evaluated in PHENIX PbSc calorimeter (lead-scintillator sandwich, ~5 x 5 cm towers, ~18 X0 depth)
- Cluster position reconstruction using logarithmic weighting of tower positions
- Energy-dependent cluster efficiency and position resolution characterized
- Published in IEEE Trans. Nucl. Sci. 47 (2000) 1982-1986

## Relevance to PPG12

- The sPHENIX EMCal clustering algorithms (RawClusterBuilderTemplate, RawClusterBuilderTopo) descend from the PHENIX clustering methodology described in this paper
- The concept of shower shape variables for photon identification (cluster width, energy ratios between NxN tower grids) traces back to this work
- PPG12's photon-ID BDT features (e11/e33, et1-4, weta_cogx, wphi_cogx) are modern implementations of the shower shape characterization pioneered here
- The cluster splitting algorithms are relevant to PPG12's use of the CLUSTERINFO_CEMC (split) vs CLUSTERINFO_CEMC_NO_SPLIT (unsplit) cluster nodes
- Already cited in the conference note (bib key: David:2000wfa)

## Note

This is a Tier 3 stub -- no arXiv PDF available. Published in IEEE TNS (2000), which does not use arXiv. May need manual enrichment.
