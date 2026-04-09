# Prompt Photon Production in Hadronic Collisions

## Overview

Prompt photons are photons produced directly in the hard scattering of partons, as opposed to secondary photons from hadron decays (primarily pi0 -> gamma gamma). They are among the cleanest probes of short-distance QCD dynamics because the photon escapes the collision without further strong-interaction fragmentation. The cross section for prompt photon production in pp collisions factorizes into parton distribution functions (PDFs), perturbatively calculable hard-scattering matrix elements, and (for the fragmentation component) photon fragmentation functions (FFs).

The PPG12 measurement of the isolated photon cross section at sqrt(s) = 200 GeV directly tests NLO pQCD predictions and constrains the gluon distribution in the proton at momentum fractions x ~ 0.04-0.35.

## Cross-Section Factorization

The inclusive prompt photon cross section in pp collisions is written as:

```
d sigma^gamma = d sigma^(direct) + d sigma^(fragmentation)

d sigma^(direct) / dpT deta = sum_{i,j} int dx1 dx2  f_i(x1, mu_F) f_j(x2, mu_F)
    * d hat{sigma}^gamma_{ij}(mu_R, mu_F, mu_ff)

d sigma^(frag) / dpT deta = sum_{i,j,k} int dx1 dx2 (dz/z^2)  f_i(x1, mu_F) f_j(x2, mu_F)
    * d hat{sigma}^k_{ij}(mu_R, mu_F, mu_ff)  D^gamma_k(z, mu_ff)
```

where f_i are parton distribution functions, d hat{sigma} are perturbative partonic cross sections, D^gamma_k are parton-to-photon fragmentation functions, and mu_R, mu_F, mu_ff are the renormalization, factorization, and fragmentation scales respectively. The variable z is the fraction of the fragmenting parton's momentum carried by the photon.

## Direct Photon Production

Direct photons couple to partons via the QED vertex at the hard-scattering level. At leading order (LO), two subprocesses contribute:

### 1. QCD Compton Scattering: qg -> q gamma

A quark from one proton scatters off a gluon from the other, producing a quark and a photon in the final state. The LO partonic cross section is:

```
d hat{sigma} / d hat{t} (qg -> q gamma) = -(alpha_em alpha_s e_q^2) / (3 hat{s}^2)
    * [ hat{s}/hat{t} + hat{t}/hat{s} ]
```

where hat{s}, hat{t} are Mandelstam variables, e_q is the quark electric charge in units of e, and the color factor is -1/3 (from the average over initial-state colors).

**This is the dominant subprocess at RHIC and LHC midrapidity**, contributing approximately 75% of the isolated photon yield at sqrt(s) = 200 GeV for pT = 8-35 GeV. The Compton process makes the prompt photon cross section directly sensitive to the gluon PDF g(x, Q^2), since one initial-state parton is always a gluon. At midrapidity, the dominant momentum fractions probed are x ~ 2 pT / sqrt(s), giving x ~ 0.08-0.35 for the PPG12 kinematic range.

### 2. Quark-Antiquark Annihilation: q qbar -> g gamma

A quark and antiquark annihilate into a gluon and a photon. The LO cross section is:

```
d hat{sigma} / d hat{t} (q qbar -> g gamma) = (8 alpha_em alpha_s e_q^2) / (9 hat{s}^2)
    * [ hat{t}/hat{u} + hat{u}/hat{t} ]
```

This process is suppressed relative to Compton scattering at RHIC energies because it requires an antiquark from the proton sea, whose PDF is much smaller than the gluon PDF at the relevant x values. Annihilation becomes more important at very high x_T (where valence quark-antiquark overlap is larger) and at ppbar colliders where the antiquark comes from the antiproton valence distribution. At the Tevatron, annihilation overtakes Compton for ET above about 120 GeV.

### NLO Corrections to Direct Production

At next-to-leading order (O(alpha_s^2)), additional subprocesses appear:

- Virtual one-loop corrections to the LO processes
- Real-emission corrections: qg -> q gamma g, q qbar -> g gamma g, qq -> qq gamma, etc.
- New channels opening at NLO: gg -> q qbar gamma (box diagram, numerically small)

The NLO K-factor (ratio of NLO to LO cross section) is typically 1.3-2.0 depending on kinematics and scale choice. NLO corrections substantially reduce the dependence on the renormalization and factorization scales, from factors of 2-3 variation at LO to approximately +/-10% at NLO when varying scales from pT/2 to 2pT.

Final-state collinear singularities arise when a quark radiates a photon collinearly. These are absorbed into the photon fragmentation functions via the MS-bar scheme, generating a factorization-scale dependence in both the direct and fragmentation components that cancels in their sum.

## Fragmentation Photons

Fragmentation photons are produced when a hard-scattered parton (quark or gluon) fragments into a photon plus hadrons. The photon carries a fraction z of the parent parton's momentum. The fragmentation contribution enters at the same formal order as the NLO direct contribution (O(alpha_s^2 alpha_em)) and cannot be neglected for a consistent NLO calculation.

At LO, the dominant fragmentation mechanism is:

```
qq -> qq  followed by  q -> gamma X  (via D^gamma_q(z, mu_ff))
```

The fragmentation functions D^gamma_q(z) and D^gamma_g(z) describe the probability for a quark or gluon to produce a photon carrying momentum fraction z. They satisfy inhomogeneous DGLAP evolution equations and are parametrized by the [BFG (Bourhis-Fontannaz-Guillet) sets](fragmentation-functions.md).

Without isolation, fragmentation photons contribute a significant fraction of the prompt photon yield, especially at moderate pT. At Tevatron kinematics (sqrt(s) = 1.96 TeV) without isolation, the fragmentation component is approximately 35% of the total at pT = 15 GeV. At sqrt(s) = 200 GeV the fraction is somewhat different due to the different parton luminosities.

## Subprocess Fractions at RHIC Kinematics

At sqrt(s) = 200 GeV and midrapidity (the PPG12 configuration), the isolated photon yield is dominated by:

| Subprocess | Approximate fraction | Notes |
|------------|---------------------|-------|
| qg -> q gamma (Compton) | ~75% | Directly sensitive to g(x) |
| q qbar -> g gamma (annihilation) | ~15% | Suppressed by sea quark PDF |
| Fragmentation (all channels) | ~10% | After isolation, see below |

These fractions are established by JETPHOX decomposition studies (Ichou & d'Enterria, arXiv:1005.4529) and are broadly consistent between sqrt(s) = 200 GeV and TeV-scale colliders at midrapidity, with Compton dominance persisting up to ET ~ 120 GeV at the Tevatron.

At forward rapidity (e.g., y = 4 at the LHC), the fragmentation fraction rises dramatically to 50-80% even after isolation, making forward photons less clean for gluon PDF extraction.

## Effect of Isolation on Production Mechanisms

[Isolation cuts](../../concepts/isolation-cuts.md) impose a maximum on the hadronic transverse energy within a cone around the photon candidate. This has a profound effect on the theoretical cross section:

### Suppression of Fragmentation

Isolation restricts the fragmentation variable z to values above a threshold:

```
z > z_c = 1 / (1 + epsilon_h)
```

where epsilon_h = E_T^max / E_T^gamma is the fractional isolation threshold. For the PPG12 parametric isolation (reco_iso_max = 0.50 + 0.043 * ET), at ET = 20 GeV the effective epsilon_h ~ 0.07, giving z_c ~ 0.93. This eliminates fragmentation photons with z < 0.93, which constitute the majority of the fragmentation contribution. The net suppression is approximately 60-80% of the fragmentation component.

After isolation with the PPG12 cone (R = 0.3, parametric cut), the residual fragmentation contribution is approximately 10-20% of the total cross section. This residual comes from photons that carry nearly all the parent parton's momentum (z close to 1), where little hadronic energy is left to violate the isolation criterion.

### Enhancement of Compton Dominance

Because isolation preferentially removes fragmentation photons while leaving direct photons largely unaffected, isolated photon measurements are even more dominated by the Compton process than inclusive ones. This enhances the sensitivity to the gluon PDF and simplifies the theoretical interpretation.

### Standard Cone vs. Smooth Cone Isolation

Two isolation prescriptions exist in the theoretical literature:

**Standard (fixed) cone** -- used experimentally by PPG12, ATLAS, CMS:
```
E_T^had < E_T^max  inside  cone R
```
This retains a nonzero fragmentation contribution for z > z_c. The NLO calculation with standard cone isolation was proven to factorize by Catani, Fontannaz, Guillet, and Pilon (hep-ph/0204023).

**Frixione smooth cone** (hep-ph/9801442):
```
E_T^had(r) < epsilon * E_T^gamma * [(1 - cos r)/(1 - cos R)]^n  for all r < R
```
This completely eliminates the fragmentation contribution (restricts it to the zero-measure set z = 1) while remaining infrared safe. It is the theoretically preferred prescription but does not exactly correspond to experimental selections.

At the PPG12 kinematics (high pT, tight isolation), the two prescriptions agree within approximately 5% (Ichou & d'Enterria, arXiv:1005.4529), so this distinction has minor practical impact.

## Physics Motivation for the PPG12 Measurement

The isolated photon cross section at sqrt(s) = 200 GeV provides:

1. **Test of NLO pQCD**: The steeply falling spectrum (approximately pT^{-7}) spans several orders of magnitude and tests the QCD hard-scattering formalism at a unique energy scale between fixed-target and collider regimes.

2. **Gluon PDF constraint**: Compton dominance at midrapidity makes the measurement directly sensitive to g(x, Q^2) at x ~ 0.04-0.35 and Q ~ 8-35 GeV. This complements DIS scaling violations and jet production data used in global PDF fits.

3. **Baseline for heavy-ion physics**: The pp cross section provides the reference for nuclear modification factor R_AA measurements in Au+Au and d+Au collisions, testing initial-state nuclear effects and final-state energy loss.

4. **x_T scaling**: When expressed as a function of x_T = 2pT/sqrt(s), prompt photon cross sections from different sqrt(s) values collapse onto a universal curve (up to logarithmic scaling violations), testing the universality of the pQCD subprocess picture.

## PPG12 Implementation

The PPG12 analysis compares the measured isolated photon cross section to two independent NLO calculations:

- **JETPHOX** (Aurenche et al., hep-ph/0602133): NLO parton-level Monte Carlo with CT14 PDFs and BFG set II fragmentation functions, run at three scales mu = 0.5 ET, ET, 2 ET. See [NLO Calculations](nlo-calculations.md).
- **Vogelsang NLO calculation** (Gordon & Vogelsang, PRD 48 (1993) 3136): Independent NLO calculation with CTEQ6M PDFs, providing a cross-check with different code and PDF set.

Both calculations include the direct and fragmentation components. The data-theory comparison is shown in `plotting/plot_final_selection.C`, which overlays JETPHOX histograms (from `NLO/MakeJetPHOXhisto.C`) and Vogelsang tabulated predictions (from `plotting/sphenix_nlo/photons_newphenix_sc*.dat`).

## Cross-References

- [NLO Calculations](nlo-calculations.md) -- JETPHOX program details, scale choices, PPG12 usage
- [Fragmentation Functions](fragmentation-functions.md) -- BFG set II parametrization and DGLAP evolution
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- Experimental isolation implementation in PPG12
- [ABCD Method](../../concepts/abcd-method.md) -- Background subtraction exploiting isolation and shower shape
- Key source papers: Owens RMP 1987 (pedagogical review), Ichou & d'Enterria arXiv:1005.4529 (subprocess fractions), Frixione hep-ph/9801442 (smooth cone isolation)
