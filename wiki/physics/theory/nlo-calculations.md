# NLO Calculations for Isolated Photon Production

## Overview

PPG12 compares its measured isolated photon cross section to two independent next-to-leading order (NLO) perturbative QCD calculations:

1. **JETPHOX** (Aurenche, Fontannaz, Guillet, Pilon, Werlen) -- the primary NLO parton-level Monte Carlo program
2. **Vogelsang NLO calculation** (Gordon & Vogelsang) -- an independent NLO code providing a cross-check

Both calculations compute the full O(alpha_s^2 alpha_em) prompt photon cross section, including direct and fragmentation components. The comparison to two independent calculations with different codes, PDFs, and fragmentation functions provides a robust test of the theoretical prediction.

## JETPHOX

### What It Calculates

JETPHOX is a general-purpose NLO partonic Monte Carlo program for prompt photon production in hadronic collisions. It computes the differential cross section:

```
d sigma^gamma / dpT deta = d sigma^(direct) + d sigma^(fragmentation)
```

at NLO accuracy, where:

- **Direct component**: Born-level processes (qg -> q gamma, q qbar -> g gamma) plus O(alpha_s) virtual and real corrections. The NLO corrections include additional radiation (qg -> q gamma g, etc.), loop diagrams, and new channels (gg -> q qbar gamma box diagram, though this is numerically small and typically switched off).

- **Fragmentation component**: A hard parton produced in a 2->2 scattering subsequently fragments into a photon. At NLO, this includes LO 2->2 parton scattering (O(alpha_s^2)) convoluted with photon fragmentation functions, plus O(alpha_s^3) real and virtual corrections to the partonic scattering.

The infrared and collinear singularities are handled via a combined phase-space slicing and subtraction method. The phase space is divided using unphysical parameters (pTm, RTh) into soft, final-state collinear, and finite-remainder regions. Singularities are regulated by dimensional continuation (n = 4 - 2 epsilon), cancelled between real and virtual contributions, and remaining collinear divergences are absorbed into PDFs and fragmentation functions in the MS-bar scheme. The final result is independent of the unphysical slicing parameters.

### Isolation Implementation

JETPHOX implements the standard fixed-cone isolation criterion:

```
E_T^had < E_T^max    inside    R = sqrt((Delta eta)^2 + (Delta phi)^2) < R_cone
```

This matches the experimental isolation used by PPG12 (and ATLAS, CMS). The isolation restricts the fragmentation variable z to z > z_c = 1/(1 + epsilon_h), where epsilon_h = E_T^max / E_T^gamma. The factorization of the isolated cross section was proven to all orders in QCD perturbation theory by Catani, Fontannaz, Guillet, and Pilon (JHEP 0205:028, 2002).

JETPHOX can also run with Frixione smooth-cone isolation, which completely eliminates the fragmentation contribution but does not correspond to the experimental procedure.

**Caveat on small cones**: The NLO calculation contains unresummed logarithms of the cone radius, log(R). For R <= 0.3 (the PPG12 cone size), these terms can become significant and the NLO calculation may not be fully reliable. The Catani et al. study found that R = 0.4 is safely in the regime where the small-cone approximation is adequate and the full-R calculation is stable, but R = 0.1 gives unphysical results (isolated cross section exceeding non-isolated). The PPG12 choice of R = 0.3 is borderline, which is an important caveat for the theory uncertainty.

### Inputs

**Parton Distribution Functions**: CT14 (Dulat et al., arXiv:1506.07443). The PPG12 final plot legend states "CT14 PDF". CT14 provides LO, NLO, and NNLO PDF sets from a global fit to 2947 data points across 33 experiments. The gluon PDF at the relevant x values (x ~ 0.04-0.35 for PPG12 kinematics) is primarily constrained by inclusive jet production at the Tevatron and LHC, and by DIS scaling violations. The gluon PDF uncertainty at Q ~ 100 GeV is approximately 3-5% in this x range.

CT14 central sets use alpha_s(M_Z) = 0.118, matching the JETPHOX default. The 56 Hessian error sets (28 eigenvector pairs at 90% C.L.) could in principle be used to compute a separate PDF uncertainty band, though PPG12 currently uses only the central set.

**Fragmentation Functions**: [BFG set II](fragmentation-functions.md) (Bourhis, Fontannaz, Guillet, Eur. Phys. J. C 2, 529, 1998). These provide D^gamma_q(z, mu_ff) and D^gamma_g(z, mu_ff) at beyond-leading-logarithm accuracy. The choice of BFG set I vs set II has negligible impact on isolated photon spectra because isolation suppresses 60-80% of the fragmentation contribution, and the two sets differ mainly in the gluon fragmentation function at low z. This was explicitly verified by Ichou & d'Enterria (arXiv:1005.4529).

**Other parameters**:
- alpha_s(M_Z) = 0.118
- N_f = 5 active quark flavors
- MS-bar renormalization and factorization scheme
- alpha_em = 1/137 (at the Thomson limit; NNLOJET uses the G_mu scheme value 1/132.232)

### Scale Choices

The three QCD scales (renormalization mu_R, initial-state factorization mu_F, fragmentation mu_ff) are set equal to a common value mu:

| Scale label | mu value | PPG12 file |
|-------------|----------|------------|
| "05" | mu = 0.5 * ET | `NLO/rootFiles/jetPHOX_05.root` |
| "10" (nominal) | mu = ET | `NLO/rootFiles/jetPHOX_10.root` |
| "20" | mu = 2.0 * ET | `NLO/rootFiles/jetPHOX_20.root` |

The scale variation from 0.5 ET to 2 ET produces approximately +/-10% variation in the cross section and is taken as the theory uncertainty band. The scale mu = pT/2 is often preferred as the central value in the literature (motivated by NLL resummation calculations and empirically providing the best fit to fixed-target and ISR data), but PPG12 uses mu = ET as the central value with the envelope from 0.5 ET to 2 ET.

### JETPHOX Output and PPG12 Processing

The JETPHOX installation is at `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/`. Output ROOT files are stored in the `pawres/` subdirectory:

- `ggdrhic_XX.root` -- direct photon component at scale XX
- `ggorhic_XX.root` -- fragmentation component at scale XX

where XX = "05", "10", or "20".

These files contain TTree `t2` with parton-level event records. The tree's UserInfo list stores the normalization vector: `norma = xsec / nb_evt / nseg`, where nseg is the number of segmentation runs (typically 100).

The macro `NLO/MakeJetPHOXhisto.C` processes these trees:

1. Reads pT bin edges from the analysis config (`pT_bins_truth`)
2. Loops over entries in both `ggdrhic` (direct) and `ggorhic` (fragmentation) trees
3. Applies rapidity cut |y| < 0.7 (matching the PPG12 eta acceptance)
4. Fills a histogram weighted by `pdf_weight[0]`
5. Scales by normalization factor and divides by bin width to produce d sigma / dpT
6. Sums direct + fragmentation into a combined histogram
7. Writes to `NLO/rootFiles/jetPHOX_XX.root`

The three scale variations are run in parallel by `NLO/run_jetphox_histos.sh`.

The final plotting macro `plotting/plot_final_selection.C` reads these histograms and constructs asymmetric error bands: the upper edge is mu = 0.5 ET (which gives a larger cross section due to larger alpha_s), and the lower edge is mu = 2 ET.

### JETPHOX Truth Isolation Studies

The macro `NLO/PlotTruthIso.C` studies the effect of the truth-level isolation cut on the JETPHOX predictions. It reads JETPHOX output for several isolation thresholds (0.5, 1, 2, 3, 4, 200 GeV) and computes the isolated cross section as a function of ET, separately for the direct and fragmentation components. This is used to validate that the JETPHOX isolation parameters match the experimental selection and to quantify the fragmentation fraction at different isolation thresholds.

## Vogelsang NLO Calculation

Gordon and Vogelsang (Phys. Rev. D 48, 3136, 1993) performed a complete NLO calculation of prompt photon production in hadronic collisions, covering both unpolarized and polarized cases. Their code is independent of JETPHOX and uses a different computational approach.

### PPG12 Usage

The Vogelsang predictions are provided as tabulated data files in `plotting/sphenix_nlo/`:

| File | Scale |
|------|-------|
| `photons_newphenix_sc1.dat` | mu = pT |
| `photons_newphenix_sc05.dat` | mu = pT/2 |
| `photons_newphenix_sc2.dat` | mu = 2 pT |

Each file contains 38 (pT, direct_yield, frag_yield, foo) entries. The plotting code reads these files, sums direct + fragmentation yields, and integrates over the PPG12 pT bins to produce bin-averaged predictions:

```cpp
myfile >> pt >> yield >> fragyield >> foo;
yield = (yield + fragyield) * factorCommon * pt;
```

The bin-averaged cross section is computed by integrating the interpolated TGraph over each bin:
```cpp
y_yield = fphoton->Integral(xmin, xmax) / (xmax - xmin);
```

Scale variations are handled identically to JETPHOX, producing an asymmetric theory band around the nominal (sc1) prediction.

### Comparison to JETPHOX

The Vogelsang calculation uses CTEQ6M PDFs (rather than CT14) and an older code base. Any differences between the JETPHOX and Vogelsang predictions arise from:

- Different PDF sets (CT14 vs CTEQ6M): approximately 2-5% effect
- Different computational methods (Monte Carlo integration vs analytic)
- Potentially different treatment of isolation (the Vogelsang files may correspond to different isolation parameters)
- Different fragmentation function implementations

The purpose of showing both on the PPG12 final plot is to demonstrate that the data-theory agreement is not an artifact of a particular code or PDF choice. The legend in `plot_final_selection.C` labels these as "NLO pQCD JETPHOX" with "CT14 PDF / BFG II FF", and "NLO pQCD by W. Vogelsang".

## Scale Dependence and Theoretical Uncertainty

The residual scale dependence at NLO is the dominant source of theoretical uncertainty. The cross section varies by approximately +/-10% when the common scale mu varies from pT/2 to 2pT, comparable for both isolated and non-isolated cases (Catani et al. found approximately 16% total variation at Tevatron kinematics).

This scale dependence arises from the truncation of the perturbative series: at NLO (O(alpha_s^2)), the uncancelled logarithms of the form alpha_s^3 ln(mu/pT) produce the residual mu dependence. Higher-order calculations (NNLO) dramatically reduce this: ATLAS (arXiv:1908.02746) showed that NNLOJET NNLO calculations reduce scale uncertainties to 0.6-5%, a factor of 2-20 improvement over NLO.

### Sources of Theoretical Uncertainty

| Source | Approximate size | Method |
|--------|-----------------|--------|
| Scale variation (mu = 0.5-2 ET) | +/-10% | Three-point variation, asymmetric band |
| PDF uncertainty (CT14 Hessian sets) | ~3-5% | Not currently propagated in PPG12 |
| alpha_s(M_Z) | ~2% | Correlated with PDFs |
| FF choice (BFG I vs II) | < 1% | Negligible after isolation |
| Isolation cone size (R ~ 0.3) | borderline | log(R) terms not resummed at NLO |
| Non-perturbative corrections | ~1% | Hadronization + underlying event |
| Missing higher orders (NNLO) | ~5-10% | Estimated from NLO/NNLO comparison at LHC |

## PHENIX Predecessor Comparison

The PHENIX measurement of direct photons at sqrt(s) = 200 GeV (arXiv:1205.5533) compared to the same Vogelsang NLO calculation with CTEQ6M PDFs and BFG II FFs. That analysis found good agreement between data and NLO predictions at scales mu = pT/2, pT, 2pT for 5.5 < pT < 25 GeV/c. PPG12 extends this comparison with dramatically improved statistics (49.6 pb^-1 vs 8 pb^-1), a wider eta acceptance (|eta| < 0.7 vs |eta| < 0.25), and both JETPHOX and Vogelsang calculations on the same plot.

## Beyond NLO

The current state of the art for isolated photon production is NNLO QCD, achieved by the NNLOJET program. ATLAS (arXiv:1908.02746) demonstrated excellent data-theory agreement at NNLO with dramatically reduced theoretical uncertainties. Key differences from NLO:

- NNLOJET computes only the direct component at NNLO; fragmentation remains NLO
- Uses hybrid-cone isolation (Frixione at small Delta R, fixed cone at large Delta R)
- Scale uncertainties reduced to 0.6-5% (vs 10-15% at NLO)
- Currently available only at sqrt(s) = 13 TeV; extension to 200 GeV would be valuable for PPG12

PPG12 compares to NLO predictions only. The NNLO comparison represents a natural future extension when NNLOJET predictions at sqrt(s) = 200 GeV become available.

## PPG12 Code References

| File | Purpose |
|------|---------|
| `NLO/MakeJetPHOXhisto.C` | Process JETPHOX output trees into pT histograms |
| `NLO/run_jetphox_histos.sh` | Run all three scale variations in parallel |
| `NLO/PlotTruthIso.C` | Study JETPHOX isolation threshold dependence |
| `NLO/rootFiles/jetPHOX_{05,10,20}.root` | Binned JETPHOX predictions at three scales |
| `plotting/sphenix_nlo/photons_newphenix_sc{05,1,2}.dat` | Vogelsang tabulated predictions |
| `plotting/plot_final_selection.C` | Final cross-section plot overlaying data + both NLO predictions |
| JETPHOX install: `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/` | JETPHOX program and raw output |

## Cross-References

- [Prompt Photon Production](prompt-photon-production.md) -- Production mechanisms and subprocess fractions
- [Fragmentation Functions](fragmentation-functions.md) -- BFG set II used in JETPHOX
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- Experimental isolation matching the JETPHOX cone criterion
- [Plotting and Systematics](../../pipeline/05-plotting-systematics.md) -- Pipeline stage that produces the data vs NLO comparison
- Key source papers: Aurenche et al. hep-ph/0602133 (JETPHOX), Catani et al. hep-ph/0204023 (NLO factorization proof), Gordon & Vogelsang PRD 48 (1993) 3136, Dulat et al. arXiv:1506.07443 (CT14)
