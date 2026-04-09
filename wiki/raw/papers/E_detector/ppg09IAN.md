---
title: "Inclusive Jet Cross-Section in sqrt(s) = 200 GeV p+p Collisions Collected in Run 2024"
authors: "H. Jiang, A. Narde, J. James, A. Brahma, J. Clement, N. Applegate, A. Hodges"
journal: "sPHENIX-Invenio (2025)"
category: E_detector
priority: 2
status: stub_needs_enrichment
tags: [sPHENIX, jet, cross-section, luminosity, IAN, internal]
ppg12_relevance: "sPHENIX internal analysis note documenting the luminosity calculation and MBD trigger cross-section used in the PPG12 cross-section normalization"
---

# PPG09 Internal Analysis Note: Inclusive Jet Cross-Section in pp at 200 GeV

## Summary

The PPG09 Internal Analysis Note (IAN) documents the sPHENIX inclusive jet cross-section measurement in pp collisions at sqrt(s) = 200 GeV using Run 2024 data. Critically for PPG12, this note contains the luminosity determination methodology: the integrated luminosity is calculated by counting minimum bias collisions for which the MBD trigger was live, corrected for prescale factors and vertex reconstruction criteria, and normalized using the measured MBD trigger cross-section. The luminosity determination is a shared infrastructure result used by multiple sPHENIX analyses.

## Key Points

- Documents the sPHENIX inclusive jet cross-section measurement (PPG09), which shares Run 2024 pp data with PPG12
- Contains the luminosity calculation methodology:
  - Prescale-corrected counting of MBD-triggered events
  - Vertex reconstruction efficiency corrections
  - MBD trigger cross-section: 25.2 +2.3/-1.7 mb in pp at sqrt(s) = 200 GeV (measured via Vernier scans)
- The luminosity result is: ~16.9 pb^-1 for the 1.5 mrad crossing angle period (run 51274-54000)
- Published as an sPHENIX internal document on the Invenio repository
- URL: https://sphenix-invenio.sdcc.bnl.gov/me/requests/ca04e089-3024-4267-aafd-81f081b4af6c

## Relevance to PPG12

- PPG12 directly uses the luminosity value from this IAN to normalize the isolated photon yield into a cross-section: d sigma / d pT = N_photon / (L_int * delta_pT * delta_eta * epsilon)
- The luminosity systematic uncertainty (~9% from MBD cross-section uncertainty) is one of the dominant systematic uncertainties in the PPG12 cross-section measurement
- The analysis note cites this as the primary luminosity reference (bib key: ppg09IAN)
- The dangling citation "MBDxsec" in systematics.tex likely refers to the same MBD cross-section measurement documented here

## Note

This is a Tier 3 stub -- no arXiv PDF available. This is an sPHENIX internal analysis note accessible only through the sPHENIX Invenio repository. May need manual enrichment with final published numbers.
