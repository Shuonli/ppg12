---
title: "PHENIX inner detectors"
authors: "M. Allen et al."
journal: "Nucl. Instrum. Meth. A 499 (2003) 549-559"
category: E_detector
priority: 2
status: stub_needs_enrichment
tags: [PHENIX, MBD, BBC, trigger, vertex, minimum-bias, Cherenkov]
ppg12_relevance: "Documents the Beam-Beam Counter (BBC, now MBD) detector design used as the minimum bias trigger and vertex detector in sPHENIX for luminosity determination and event selection in PPG12"
---

# PHENIX Inner Detectors

## Summary

This paper, part of the NIM A 499 RHIC special volume, describes the PHENIX inner detector systems, with particular focus on the Beam-Beam Counters (BBC). The BBC consists of two arrays of 64 quartz Cherenkov radiators coupled to mesh-dynode photomultiplier tubes, positioned at forward and backward rapidity (3.0 < |eta| < 3.9) on either side of the interaction point. The BBC provides the minimum bias trigger, collision vertex determination (via time-of-flight difference between north and south arms), and centrality measurement for PHENIX. The detector was later inherited by sPHENIX and renamed the Minimum Bias Detector (MBD), where it continues to serve the same functions.

## Key Points

- BBC/MBD design: 2 x 64 photomultiplier tubes with quartz Cherenkov radiators
- Positioned at +/- 144 cm from the interaction point, covering 3.0 < |eta| < 3.9
- Timing resolution: ~50-100 ps per PMT, enabling vertex determination to ~2 cm via (t_north - t_south) * c / 2
- Provides Level-1 minimum bias trigger based on coincidence of hits in north and south arms
- Charge sum provides centrality determination in heavy-ion collisions
- Timing spread (sigma_t) across PMTs used as a pileup metric in sPHENIX (MbdPileupHelper.h)
- Trigger cross-section in pp at 200 GeV: 25.2 +2.3/-1.7 mb (measured by sPHENIX via Vernier scans)

## Relevance to PPG12

- The MBD (ex-BBC) is the primary trigger detector for PPG12 minimum bias data collection
- MBD timing is used for vertex reconstruction, with |vtx_z| < 10 cm as a standard event selection cut
- MBD timing spread (avg sigma_t from north and south arms) is used as a pileup metric:
  - `mbd_avg_sigma_max` cut in analysis configs rejects high-pileup events
  - MbdPileupHelper.h computes timing metrics from 128 PMT channels
  - Events with multiple collisions per crossing show broader timing distributions
- The MBD trigger cross-section determines the integrated luminosity via L = N_MB / sigma_MBD
- Already cited in the conference note (bib key: phenix_mbd_nim)

## Note

This is a Tier 3 stub -- no arXiv PDF available. Published in the RHIC special volume of NIM A (2003), which predates routine arXiv posting of instrumentation papers. May need manual enrichment.
