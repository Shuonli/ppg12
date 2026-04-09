# Photon Fragmentation Functions: BFG Set II

## Overview

Photon fragmentation functions D^gamma_a(z, M^2) describe the probability for a parton of flavor a (quark or gluon) to produce a photon carrying a fraction z of the parton's momentum, at a factorization scale M. They enter the fragmentation component of the prompt photon cross section and are required for any NLO calculation that uses standard (fixed-cone) isolation.

PPG12 uses the **BFG set II** fragmentation functions (Bourhis, Fontannaz, Guillet, Eur. Phys. J. C 2, 529, 1998; arXiv:hep-ph/9704447) in its JETPHOX NLO predictions. The final plot legend explicitly reads "CT14 PDF / BFG II FF".

## Structure of the Photon Fragmentation Function

The full fragmentation function decomposes into two physically distinct pieces:

```
D^gamma(z, M^2) = D^{gamma,AN}(z, M^2) + D^{gamma,NP}(z, M^2)
```

### Anomalous (Perturbative) Component: D^{gamma,AN}

The anomalous component is fully calculable in perturbative QCD. It arises from the pointlike coupling of the photon to quarks (the QED vertex gamma -> q qbar), which generates an inhomogeneous driving term in the DGLAP evolution equations. At leading-logarithm (LL) accuracy, D^{gamma,AN} grows as:

```
D^{gamma,AN} ~ (4 pi alpha_em / alpha_s) * ln(M^2 / Lambda_QCD^2)
```

This is unique to the photon: unlike hadron fragmentation functions (which are entirely non-perturbative at the initial scale), the photon FF has a perturbatively calculable piece that grows logarithmically with scale. The anomalous component dominates at large M^2 and moderate z.

### Non-Perturbative Component: D^{gamma,NP}

The non-perturbative component encodes the long-distance physics at the initial scale M_0^2 ~ m_rho^2 ~ 0.5 GeV^2. It is modeled using Vector Meson Dominance (VDM) -- the idea that the photon fluctuates into vector mesons (rho, omega, phi) that then fragment like hadrons:

```
gamma = g * sqrt(2) * [ rho + omega/3 - sqrt(2)/3 * phi ]
```

The non-perturbative fragmentation is extracted from rho-meson fragmentation data measured at LEP (ALEPH) and at lower-energy e+e- experiments (HRS):

```
D^{gamma,VDM}_q = g^2 * [ (4/9) D^{(u ubar)}_q + (1/9) D^{(d dbar)}_q + (1/9) D^{(s sbar)}_q ]
```

where g^2 ~ alpha_em and D^{(q qbar)}_a are the rho-meson fragmentation functions. The non-perturbative component satisfies the standard (homogeneous) DGLAP evolution for M^2 > M_0^2.

## DGLAP Evolution Equations

The fragmentation functions obey inhomogeneous DGLAP equations. In the singlet sector:

```
M^2  dD^gamma_q / dM^2 = C_s * K^gamma_q  +  P_qq (x) D^gamma_q  +  P_gq (x) D^gamma_g

M^2  dD^gamma_g / dM^2 = C_s * K^gamma_g  +  P_qg (x) D^gamma_q  +  P_gg (x) D^gamma_g
```

where:
- P_ij are the standard timelike splitting functions (same as for hadron FFs)
- (x) denotes Mellin convolution
- K^gamma_a are the inhomogeneous kernels specific to the photon
- C_s is a color factor

The inhomogeneous kernels have the perturbative expansion:

```
K^gamma_a(z, M^2) = (alpha_em / 2 pi) * [ K^(0)_{gamma a}(z) + (alpha_s / 2 pi) K^(1)_{gamma a}(z) + ... ]
```

The leading kernel for quarks is the QED splitting function:

```
K^(0)_{gamma q}(z) = e_q^2 * P_{gamma q}(z) = e_q^2 * [1 + (1-z)^2] / z
```

The gluon couples to the photon only through a quark loop, so K^gamma_g starts at O(alpha_s):

```
K^gamma_g = O(alpha_em * alpha_s)
```

This means the gluon-to-photon fragmentation is suppressed by one power of alpha_s relative to the quark channel.

### Beyond Leading Logarithm (BLL) Accuracy

The BFG calculation goes beyond LL by including the O(alpha_s) corrections to the inhomogeneous kernels (K^(1) terms). This improves the accuracy of the anomalous component, particularly at small z where the LL approximation breaks down. At BLL, the gluon anomalous fragmentation function D^{gamma,AN}_g can become negative at very small z (z < 0.01) due to kernel singularities, but this region is not probed experimentally.

## Factorization Scheme Dependence

The non-perturbative initial condition in the MS-bar scheme includes a perturbative correction beyond the VDM model:

```
D^{gamma,NP}_{ns,i}(z, M_0^2) = D^{gamma,VDM}_{ns,i}(z, M_0^2) - C_{ns,i} * D^{gamma,MS}(z)
```

where the MS-bar correction is:

```
D^{gamma,MS}(z) = (alpha_em / 2 pi) * { [1 + (1-z)^2] / z * [ln(1-z) + ln(z)] + z }
```

This subtraction ensures factorization scheme invariance of the physical cross section: the scheme-dependent pieces in the fragmentation function cancel against scheme-dependent pieces in the partonic coefficient functions. The physical observable d sigma^gamma is independent of the choice of scheme, but the decomposition into "direct" and "fragmentation" components is scheme-dependent.

## BFG Set I vs Set II

Two parameterizations of the VDM initial condition are provided, differing primarily in the gluon-to-photon fragmentation:

### VDM Parametrization at Q_0^2 = 2 GeV^2

The rho-meson fragmentation functions are parametrized as D = N * x^alpha * (1-x)^beta:

| Component | Set I (N, alpha, beta) | Set II (N, alpha, beta) |
|-----------|----------------------|------------------------|
| Valence | (0.785, -0.5, 1.499) | (1.140, -0.2, 1.693) |
| Sea | (0.111, -1, 2.912) | (0.100, -0.3, 3.437) |
| Gluon | (0.108, -1, 3) | (2.550, -0.3, 3) |

The key difference is the gluon normalization: Set II ("large gluon") has N_gluon = 2.55 vs Set I ("small gluon") with N_gluon = 0.108, a factor of ~24 difference. This is because Set II constrains the ratio D^rho_g / D^rho_u to be comparable to D^pi0_g / D^pi0_u, while Set I leaves the gluon less constrained.

- **Set I**: chi^2/dof = 1.33 (fit to ALEPH and HRS data)
- **Set II**: chi^2/dof = 1.22 (slightly better fit)

### Impact on Isolated Photon Predictions

Despite the large difference in gluon FFs, the isolated photon cross section is virtually identical for BFG I and BFG II. This was explicitly demonstrated by Ichou & d'Enterria (arXiv:1005.4529). The reason is twofold:

1. **Isolation suppresses fragmentation by 60-80%**, reducing the overall sensitivity to FFs.
2. The gluon FF enters at O(alpha_s) relative to the quark FF, making D^gamma_g subdominant even for the fragmentation component.
3. The gluon FF difference is largest at low z, which is precisely the region excluded by the isolation cut (z > z_c ~ 0.93 for tight isolation).

PPG12 uses Set II as the default, following the standard choice in the JETPHOX literature.

## Heavy Quark Treatment

Heavy quark (charm, bottom) fragmentation functions are activated at their mass thresholds:

```
D^gamma_c(mu^2 < m_c^2) = 0    (m_c = 1.5 GeV)
D^gamma_b(mu^2 < m_b^2) = 0    (m_b = 4.5 GeV)
```

At the mass threshold, the MS-bar input is:

```
D^{gamma,MS}_Q(z, m_Q^2) = -(alpha_em / 2 pi) * e_Q^2 * [1 + (1-z)^2] / z * (2 ln z + 1)
```

The QCD scale is Lambda^(4)_QCD = 230 MeV. Above the threshold, heavy quark FFs evolve according to the standard DGLAP equations with N_f + 1 active flavors.

## Role in the NLO Cross Section

In the JETPHOX NLO calculation, the fragmentation component takes the form:

```
d sigma^(frag) / dpT deta = sum_{i,j,k} int dx1 dx2 (dz/z^2) f_i(x1, mu_F) f_j(x2, mu_F)
    * d hat{sigma}^k_{ij}(mu_R, mu_F, mu_ff) * D^gamma_k(z, mu_ff)
```

The fragmentation scale mu_ff controls the evolution of D^gamma_k. In the PPG12 JETPHOX runs, mu_ff = mu_R = mu_F = mu (common scale), varied from 0.5 ET to 2 ET.

The mu_ff dependence of the fragmentation function is:

```
d ln D^gamma_q / d ln mu_ff^2 ~ (alpha_em / 2 pi) * e_q^2 * P_{gamma q} + (alpha_s / 2 pi) * [P_qq (x) D^gamma_q + P_gq (x) D^gamma_g] + ...
```

This scale dependence partially cancels against the mu_ff dependence of the partonic coefficient function d hat{sigma}^k, so the physical cross section has reduced sensitivity to mu_ff. The residual dependence is part of the +/-10% scale uncertainty band.

### Fragmentation Fraction at RHIC Kinematics

At sqrt(s) = 200 GeV with standard cone isolation (R = 0.3, parametric ET-dependent cut), the fragmentation contribution to the total isolated photon cross section is approximately:

- ~10-20% at pT = 10 GeV (moderate isolation suppression)
- ~5-10% at pT = 30 GeV (stronger effective isolation at higher pT)

This can be studied directly with the `NLO/PlotTruthIso.C` macro, which reads JETPHOX output for different isolation thresholds and decomposes the cross section into direct (`ggdrhic`) and fragmentation (`ggorhic`) components.

## Suppression Under Isolation

The mechanism by which isolation suppresses the fragmentation contribution is central to understanding why BFG set choice matters so little for PPG12.

With standard cone isolation requiring E_T^had < E_T^max inside cone R, a fragmenting parton with the photon carrying fraction z must satisfy:

```
(1 - z) * E_T^parton < E_T^max
=> z > z_c = E_T^gamma / (E_T^gamma + E_T^max) = 1 / (1 + epsilon_h)
```

For the PPG12 parametric isolation at ET = 20 GeV: E_T^max ~ 0.50 + 0.043 * 20 = 1.36 GeV, giving epsilon_h ~ 0.068 and z_c ~ 0.936. The fragmentation function is then integrated only over z in [0.936, 1], excluding the bulk of the distribution.

In this restricted z range, D^gamma_q(z) is dominated by the perturbative (anomalous) component, which at z -> 1 behaves as:

```
D^{gamma,AN}_q(z -> 1) ~ (alpha_em / 2 pi) * e_q^2 * ln(1 / (1-z))
```

This logarithmic divergence is integrable and is cancelled in the physical cross section by the direct coefficient function (a consequence of the KLN theorem). The non-perturbative (VDM) component, which differs between BFG I and II, is concentrated at moderate z and is almost entirely excluded by the isolation cut.

## Comparison: Frixione Isolation and Fragmentation

Under Frixione smooth-cone isolation, the fragmentation contribution is restricted to the zero-measure set z = 1 and vanishes identically. This eliminates any dependence on fragmentation functions. The Frixione prescription is therefore the theoretically cleanest approach, but it does not match the experimental isolation criterion.

The fact that standard-cone and Frixione-cone isolated photon predictions agree within approximately 5% at the PPG12 kinematics (high pT, tight isolation, midrapidity) confirms that the residual fragmentation contribution is small and well-controlled.

## Cross-References

- [Prompt Photon Production](prompt-photon-production.md) -- Production mechanisms where FFs enter
- [NLO Calculations](nlo-calculations.md) -- JETPHOX program that uses BFG II
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- Experimental isolation that suppresses fragmentation
- Key source papers: Bourhis, Fontannaz, Guillet hep-ph/9704447 (BFG paper), Catani et al. hep-ph/0204023 (factorization proof with standard isolation), Ichou & d'Enterria arXiv:1005.4529 (BFG I vs II comparison)
