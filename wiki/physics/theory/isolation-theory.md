# Theoretical Foundations of Photon Isolation

## Motivation: Why Isolate Photons?

Prompt photon production in hadronic collisions proceeds through two mechanisms:

1. **Direct production** -- the photon participates in the hard scatter (e.g. qg -> gamma q Compton, qqbar -> gamma g annihilation). The photon emerges isolated from hadronic activity.
2. **Fragmentation production** -- a hard parton fragments into a photon carrying momentum fraction z, accompanied by collinear hadrons. The photon is embedded in a jet.

At NLO, the inclusive prompt photon cross section decomposes as (Catani et al., hep-ph/0204023, Eq. 1):

```
d sigma = d sigma^dir + d sigma^frag
        = sum_{a,b} int dx_a dx_b F_a(x_a) F_b(x_b)
          x [ d hat{sigma}^gamma_ab  +  sum_c int (dz/z^2) d hat{sigma}^c_ab D^gamma_c(z; mu_ff) ]
```

where F_a are parton distribution functions, D^gamma_c are parton-to-photon fragmentation functions, and the sum runs over initial-state parton flavors a,b and fragmenting parton flavor c. The fragmentation component is non-perturbative and introduces model dependence through the fragmentation functions (BFG-I, BFG-II sets).

Isolation cuts suppress fragmentation by requiring low hadronic activity around the photon. This simplifies the theoretical calculation (fewer non-perturbative inputs) and experimentally suppresses the dominant pi0/eta decay background. The theoretical challenge is imposing isolation without spoiling the infrared (IR) safety of the perturbative calculation.

## Standard (Fixed) Cone Isolation

The experimentally standard criterion requires that the total hadronic transverse energy inside a cone of radius R around the photon satisfy:

```
E_T^had(R) = sum_i E_T^i theta(R - R_{i gamma}) <= E_T^max
```

where R_{i gamma} = sqrt((eta_i - eta_gamma)^2 + (phi_i - phi_gamma)^2) and E_T^max is a fixed threshold (or proportional to E_T^gamma via epsilon_h: E_T^max = epsilon_h * E_T^gamma).

**Factorization holds.** Catani, Fontannaz, Guillet, and Pilon (hep-ph/0204023) proved to all orders that the inclusive factorized structure carries over to the isolated case. The isolation constraint modifies the fragmentation integration range: only fragmentation at z > z_c = 1/(1 + epsilon_h) contributes. The isolated cross section factorizes as:

```
sigma^iso(p_gamma; z_c, R) = sum_a int_{z_c}^1 (dz/z) hat{sigma}^{a,iso} D^gamma_a(z; mu_ff)
                              + hat{sigma}^{gamma,iso}
```

The key property is that the isolation function satisfies four conditions ensuring IR and collinear safety: (i) it is a valid measurement function on the final state, (ii) soft partons do not change isolation status, (iii) initial-state collinear splittings factorize, and (iv) final-state collinear splittings inside the cone preserve the fragmentation structure with modified z range.

**The problem.** Standard cone isolation retains a residual fragmentation contribution from z > z_c. For epsilon_h = 0.1, this means z > 0.91 -- the photon must carry at least 91% of the fragmenting parton's momentum, but the fragmentation function in this region is non-zero and poorly constrained. At Tevatron kinematics (sqrt(s) = 1.96 TeV), Catani et al. find fragmentation is roughly 12% of the total NLO isolated cross section at R = 0.7, epsilon_h = 0.133.

**Small-cone behavior.** The NLO calculation introduces logarithms of the cone radius: terms proportional to log(R) appear in the direct component when a parton is radiated just outside the cone. For R >= 0.3-0.4, the small-R approximation is adequate (overestimates fragmentation by ~8% at R = 1.0). Below R ~ 0.1 the NLO calculation becomes unreliable -- the isolated cross section can exceed the non-isolated one, signaling the need for log(R) resummation. The Catani et al. numerical results at Tevatron kinematics (Table 1):

| R    | Direct NLO (pb/GeV) | Frag NLO | Total NLO |
|------|---------------------|----------|-----------|
| 1.0  | 3318                | 447      | 3765      |
| 0.7  | 3603                | 495      | 4098      |
| 0.4  | 3969                | 556      | 4525      |
| 0.1  | 4758                | 679      | 5431*     |

*R = 0.1 gives unphysical result (isolated > non-isolated = 5218), indicating breakdown.

**Scale dependence.** The isolated cross section varies by ~16% when the common scale mu varies from p_T/2 to 2p_T, comparable to the non-isolated case. Isolation does not induce additional scale sensitivity.

## Frixione Smooth Cone Isolation

Frixione (hep-ph/9801442) introduced a smooth isolation criterion that completely eliminates the fragmentation contribution while preserving IR safety to all orders.

### The criterion

For every sub-cone of radius delta <= delta_0 around the photon, the total hadronic transverse energy must satisfy:

```
sum_i E_T^i theta(delta - R_{i gamma}) <= X(delta)    for all delta <= delta_0
```

with the threshold function:

```
X(delta) = E_T^gamma epsilon_gamma [ (1 - cos delta) / (1 - cos delta_0) ]^n
```

The canonical choice is epsilon_gamma = 1, n = 1. The defining property is X(delta) -> 0 as delta -> 0: the allowed hadronic energy vanishes smoothly at zero angular separation from the photon.

### Why it works: simultaneous IR safety and fragmentation elimination

**Soft gluons** (E -> 0): The isolation condition X(delta) provides no constraint since E < X(delta) for any delta when E is sufficiently small. Soft gluon emission is unconstrained, and the usual QCD IR cancellation between real and virtual contributions proceeds normally.

**Collinear quarks** (q nearly parallel to gamma, delta -> 0): The smooth criterion forces E_quark -> 0 proportionally to delta^{2n}. This energy damping cancels the 1/(1-y) collinear divergence. The finite contribution from inside the cone scales as:

```
sigma_cone ~ E_gamma^2 epsilon_gamma^2 / (4n)
```

which is finite for n >= 1/2.

**Collinear gluons** (g nearly parallel to gamma): The finite gluon contribution inside the cone is:

```
sigma_cone ~ (1 - cos delta_0) [ log(E_gamma epsilon_gamma) - n ]
```

For standard cone isolation (n = 0, epsilon_gamma = epsilon_c << 1), this gives a large negative correction proportional to log(epsilon_c). With n = 1 and epsilon_gamma = 1, no such large logarithm appears.

**Fragmentation suppression.** Because X(delta) -> 0 as delta -> 0, a collinear quark fragmenting into a photon at z < 1 must carry vanishing energy -- the fragmentation contribution is restricted to the zero-measure point z = 1, where it vanishes. The isolated cross section depends only on the direct production mechanism.

### Comparison: standard vs smooth cone

| Property | Standard cone | Frixione smooth cone |
|----------|--------------|---------------------|
| IR safe | Yes (if E_T^max > 0) | Yes |
| Fragmentation eliminated | No (z > z_c survives) | Yes (restricted to z = 1) |
| Large logarithmic corrections | log(epsilon_c) when epsilon_c << 1 | Absent (moderate corrections) |
| Experimental implementation | Straightforward | Requires testing all sub-cones (impractical at detector level, approximated) |
| JETPHOX support | Yes | Yes |

### Experimental practicalities

The Frixione criterion requires testing the isolation condition at every angular scale from 0 to delta_0 simultaneously. This is experimentally impractical at the calorimeter level -- detector granularity and noise set a minimum resolvable angular scale, below which X(delta) falls below the detector energy threshold. In practice, the contribution from z > z_th ~ 1 - E_th/E_gamma is negligibly small, so the distinction between smooth and standard cone isolation is quantitatively minor at high E_T.

NNLOJET (the NNLO direct-photon code) uses a "hybrid-cone" isolation: Frixione smooth cone at small delta (< 0.1) and standard fixed cone at the full radius R = 0.4, combining theoretical cleanliness with practical matching to experimental cuts.

## Quantitative Comparison: Standard vs Smooth at Collider Kinematics

Ichou and d'Enterria (1005.4529) performed a systematic JETPHOX comparison of the two isolation prescriptions:

- **Tevatron (sqrt(s) = 1.96 TeV, midrapidity):** Smooth-cone isolation gives ~5% lower cross sections than standard cone at low E_T, converging at high E_T where fragmentation is suppressed kinematically.
- **LHC (sqrt(s) = 14 TeV, midrapidity):** Similar ~5% difference at moderate E_T.
- **LHC forward rapidity (y = 4):** The difference reaches 10-15% at low E_T because fragmentation contributes 50-80% of inclusive prompt photon production in this regime, and standard cone isolation is less effective at suppressing it.

The smallness of the difference at midrapidity and high E_T validates the use of standard cone isolation for precision measurements like PPG12.

### Subprocess dominance after isolation

The isolation cut reshapes the subprocess composition. Ichou and d'Enterria show (Fig. 3 of 1005.4529):

- **Compton (qg -> gamma q):** Dominates at midrapidity for all accessible E_T at both RHIC (sqrt(s) = 200 GeV) and LHC energies. This makes isolated photon measurements directly sensitive to the gluon PDF g(x, Q^2).
- **Annihilation (qqbar -> gamma g):** Sub-dominant at midrapidity; contributes more at high xT where valence quarks dominate.
- **Fragmentation:** Suppressed from ~40% (inclusive) to ~10-20% (isolated) at midrapidity. At LHC forward rapidity, fragmentation remains 50-80% even after isolation.

At RHIC sqrt(s) = 200 GeV and midrapidity, the dominant x values probed are x ~ 2 E_T / sqrt(s) ~ 0.08-0.35 for PPG12's pT range of 8-35 GeV.

## Connection to PPG12

### PPG12's isolation choice

PPG12 uses standard cone isolation with R = 0.3 and a parametric (ET-dependent) threshold:

```
reco_iso_max = b + s * ET
```

with nominal parameters b = 0.502 GeV, s = 0.0433 (from config_bdt_nom.yaml). At ET = 20 GeV, this gives E_T^max ~ 1.37 GeV, corresponding to epsilon_h ~ 0.07 -- tighter than the typical ATLAS/CMS choice of epsilon_h ~ 0.1-0.15.

The ET-dependent threshold accounts for the growth of underlying event and pileup energy with harder interactions, maintaining a roughly constant signal efficiency as a function of ET. This is conceptually similar to the ATLAS parametric isolation (E_T^iso < 4.8 + 0.0042*ET GeV) though with different parameters optimized for sPHENIX kinematics.

### Theoretical implications

1. **Residual fragmentation is small.** At PPG12 kinematics (tight isolation, midrapidity, moderate ET), the fragmentation contribution is roughly 10-15% of the total NLO cross section. This is included in the JETPHOX NLO prediction that PPG12 compares against.

2. **Standard vs smooth difference is ~5%.** Within the PPG12 ET range and the current experimental uncertainties, the choice of standard vs Frixione isolation has negligible impact. The JETPHOX calculation uses standard cone isolation matching the experimental procedure.

3. **R = 0.3 is borderline for NLO reliability.** Catani et al. note that their NLO calculation becomes unreliable for R <= 0.3. At R = 0.3, log(R) corrections are still moderate but near the edge of applicability. PPG12's use of R = 0.3 is a compromise between background rejection (favoring smaller R) and NLO reliability (favoring larger R). Systematic assessment of the cone-size dependence is warranted.

4. **BFG fragmentation function choice is irrelevant.** Ichou and d'Enterria demonstrate that BFG-I and BFG-II give virtually identical isolated photon spectra because isolation suppresses the fragmentation contribution where the two sets differ most (low z). PPG12's use of BFG-II introduces no significant theoretical uncertainty.

5. **Scale variation prescription.** The ~10-16% cross-section variation from mu = E_T/2 to 2E_T motivates PPG12's standard practice of running JETPHOX at three scales to define the NLO theory uncertainty band.

### Truth isolation in PPG12

At the generator (truth) level, PPG12 applies a truth isolation cut of iso_ET_truth < 4 GeV in a cone R = 0.3 around the truth photon. This defines the fiducial cross section to which the measurement is corrected. The 4 GeV threshold is loose enough to retain essentially all direct photons while removing most fragmentation photons at z << 1.

## Key References

| Reference | Content |
|-----------|---------|
| Frixione, hep-ph/9801442 | Smooth cone isolation: criterion, IR safety proof, fragmentation elimination |
| Catani et al., hep-ph/0204023 | NLO standard cone isolation: factorization proof, full R dependence, JETPHOX precursor |
| Ichou & d'Enterria, 1005.4529 | Standard vs smooth quantification (~5%), subprocess decomposition, PDF sensitivity |
| Aurenche et al., hep-ph/0602133 | JETPHOX validation paper |
| Bourhis, Fontannaz, Guillet, hep-ph/9704447 | BFG fragmentation function sets |

## See Also

- [Isolation Cuts](../../concepts/isolation-cuts.md) -- PPG12-specific implementation: cone sizes, tower vs topo, parametric thresholds, MC corrections
- [Isolation Methods Across Experiments](../techniques/isolation-methods.md) -- Experimental implementations at ATLAS, CMS, ALICE, PHENIX, sPHENIX
- [ABCD Sideband Method](../techniques/abcd-sideband-method.md) -- Background subtraction using isolation as one axis
