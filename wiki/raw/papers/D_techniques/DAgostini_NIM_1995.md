---
title: "A multidimensional unfolding method based on Bayes' theorem"
authors: "G. D'Agostini"
journal: "Nucl. Instrum. Meth. A 362 (1995) 487-498"
category: D_techniques
priority: 3
status: stub_needs_enrichment
tags: [unfolding, bayesian, iterative, statistical-methods]
ppg12_relevance: "Foundational paper for the iterative Bayesian unfolding method used in PPG12 to correct the measured photon pT spectrum for detector effects"
---

# A Multidimensional Unfolding Method Based on Bayes' Theorem

## Summary

This 1995 paper by G. D'Agostini introduced the iterative Bayesian unfolding method for correcting measured distributions for detector smearing and acceptance effects. The method applies Bayes' theorem to invert the detector response matrix, using the measured distribution as a prior that is iteratively updated. It provides a regularized alternative to direct matrix inversion, which is numerically unstable when the response matrix has large off-diagonal elements. The approach became one of the two standard unfolding methods in high-energy physics (alongside SVD unfolding).

## Key Points

- Introduces the application of Bayes' theorem to the unfolding (deconvolution) problem in particle physics
- The measured spectrum is treated as a conditional probability: P(effect|cause) is inverted to obtain P(cause|effect) via Bayes' theorem
- Iterative procedure: the unfolded result from iteration N becomes the prior for iteration N+1, converging toward the true distribution
- Naturally handles multidimensional problems (multiple observables simultaneously)
- Provides built-in regularization through the choice of number of iterations -- fewer iterations = more regularization, more iterations = closer to unregularized inversion
- The method preserves the total number of events and handles bin-to-bin migration in a statistically principled way
- Later improved and extended in D'Agostini's 2010 paper (arXiv:1010.0632), which is the version most commonly implemented in modern software (RooUnfold)

## Relevance to PPG12

- PPG12 uses iterative Bayesian unfolding (via RooUnfold) to correct the reconstructed isolated photon pT spectrum for EMCal energy resolution, bin migration, and acceptance effects
- The response matrix is built from PYTHIA + GEANT4 simulation of sPHENIX, mapping truth pT to reconstructed pT
- The number of iterations is a tunable parameter that affects the bias-variance tradeoff; PPG12 optimizes this as part of the unfolding systematic uncertainty
- This 1995 paper should be cited alongside the 2010 improved version (already cited) as the original foundational reference

## Note

This is a Tier 3 stub -- no arXiv PDF available. Content is based on metadata and domain knowledge. May need manual enrichment.
