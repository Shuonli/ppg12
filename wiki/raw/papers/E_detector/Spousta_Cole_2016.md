---
title: "Interpreting single jet measurements in Pb+Pb collisions at the LHC"
authors: "Martin Spousta, Brian Cole"
journal: "Eur. Phys. J. C 76 (2016) 50"
category: E_detector
priority: 2
status: stub_needs_enrichment
tags: [modified-power-law, jet-spectrum, fitting, energy-loss, heavy-ion]
ppg12_relevance: "Provides the modified power law functional form used in PPG12 to fit the MC photon pT spectrum for smooth interpolation and ratio construction"
---

# Interpreting Single Jet Measurements in Pb+Pb Collisions at the LHC

## Summary

Spousta and Cole introduced a modified power-law parameterization for fitting steeply falling pT spectra of jets (and by extension, photons) in hadronic collisions. The paper's primary focus is on interpreting jet suppression measurements in Pb+Pb collisions at the LHC, but its key technical contribution -- the modified power-law functional form -- has been widely adopted for fitting any steeply falling particle spectrum where a simple power law is insufficient. The function accounts for the curvature of the spectrum in log-log space that arises from the running of alpha_s and PDF evolution, providing a smooth analytic description across a wide pT range with only a few free parameters.

## Key Points

- Introduces a modified power-law function for fitting particle pT spectra: d sigma / d pT = A * pT^(-n) * (1 + pT/p0)^(-m), or similar multi-parameter forms that capture the curvature beyond a simple power law
- Motivation: a pure power-law fit (d sigma / d pT ~ pT^-n) fails across wide pT ranges because the effective exponent n varies with pT due to PDF evolution and alpha_s running
- The modified form provides a smooth, well-behaved fit that can be used for:
  - Interpolation between measured data points
  - Ratio construction (data/fit) for residual analysis
  - Bin-width corrections and spectrum unfolding cross-checks
- Applied to LHC jet spectra to extract jet energy loss parameters in QGP, but the fitting methodology is general
- Published in Eur. Phys. J. C 76 (2016) 50 (also listed as EPJC 76 (2016) 2 in some indexing)

## Relevance to PPG12

- PPG12 uses the Spousta-Cole modified power-law function to fit the combined PYTHIA photon MC spectrum (photon5 + photon10 + photon20 samples, cross-section weighted) for smooth interpolation
- The fit is used in the selection chapter of the analysis note to show the MC/fit ratio, validating the smooth combination of MC samples across pT bins
- Already cited in the analysis note (bib key: Spousta_2016)

## Note

This is a Tier 3 stub -- no arXiv PDF available in the wiki collection. The paper may have an arXiv preprint (likely 1504.01945 or similar) but was not included in the PDF download batch. May need manual enrichment.
