---
title: "sPHENIX Technical Design Report, PD-2/3 Release"
authors: "The sPHENIX Collaboration"
journal: "sPHENIX internal document (2019)"
category: E_detector
priority: 3
status: stub_needs_enrichment
tags: [sPHENIX, detector, TDR, EMCal, HCal, magnet, design]
ppg12_relevance: "Definitive reference for sPHENIX detector design parameters including EMCal energy resolution, tower segmentation, solenoid field, and calorimeter geometry used throughout the PPG12 analysis"
---

# sPHENIX Technical Design Report

## Summary

The sPHENIX Technical Design Report (TDR) is the comprehensive design document for the sPHENIX detector at RHIC, released in May 2019 as the PD-2/3 milestone document. It specifies the full detector design including the tungsten-scintillating fiber (W/SciFi) electromagnetic calorimeter, the inner and outer hadronic calorimeters, the recycled BaBar superconducting solenoid providing a 1.5 T field, the compact TPC tracker, and the minimum bias detector (MBD) trigger system. The TDR provides the engineering specifications, expected performance figures, and physics projections that define the baseline detector configuration used for all sPHENIX analyses including PPG12.

## Key Points

- **EMCal**: Tungsten powder / scintillating fiber (W/SciFi) sampling calorimeter with projective tower geometry
  - Energy resolution: ~12%/sqrt(E) + 2% constant term (design specification)
  - Tower segmentation: 0.024 x 0.024 in (eta, phi), 24 x 256 = 6144 towers per sector
  - ~18 radiation lengths deep, covering |eta| < 1.1
- **Inner HCal (IHCal)**: Aluminum absorber + scintillating tiles, inside solenoid bore
- **Outer HCal (OHCal)**: Steel absorber + scintillating tiles, outside solenoid, serving as magnetic flux return
  - Combined HCal resolution: ~100%/sqrt(E)
- **Magnet**: Recycled BaBar superconducting solenoid from SLAC, providing 1.5 T uniform field
- **Tracking**: Compact TPC with GEM readout for charged particle tracking, plus MVTX and INTT silicon detectors
- **MBD (formerly BBC)**: Minimum bias trigger and vertex detector using Cherenkov radiators + PMTs
- **Acceptance**: Full azimuth, |eta| < 1.1 for calorimeters, |eta| < 1.0 for tracking
- Available at: https://indico.bnl.gov/event/7081/attachments/25527/38284/sphenix%20tdr%2020190513.pdf

## Relevance to PPG12

- Provides the baseline detector specifications used in PPG12: EMCal resolution and segmentation determine shower shape variables, energy scale, and isolation cone performance
- The EMCal tower size (0.024 x 0.024) defines the granularity of the cluster shape variables (e11/e33, et1-4, weta, wphi) used in the BDT photon identification
- The solenoid field (1.5 T) determines the magnetic bending of charged particles that contribute to the isolation energy
- GEANT4 simulations of the detector use the TDR geometry as the baseline configuration
- Already cited in both the conference note and analysis note (bib key: Sphenix:TDR)

## Note

This is a Tier 3 stub -- no arXiv PDF available. The TDR is an internal sPHENIX document hosted on the BNL Indico server. May need manual enrichment.
