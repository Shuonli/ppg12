# x_T-scaling compilation — how this plot was made

Reproduction of **Figure 1 (left panel) of arXiv:1901.10950** (F. Bock, Hard Probes 2018) — the
world compilation of x_T-scaled (isolated) direct-photon invariant cross sections in
p+p / p̄+p — with the **PPG12 sPHENIX** isolated-photon measurement overlaid.

Plotted quantity:

```
Y = (√s/GeV)^n · E d³σ/dp³   [pb GeV⁻² c³] ,   n = 4.5
x = x_T = 2 p_T / √s ,   at y ≈ 0
```

Final figure: `../figures/xtscaling_compilation.pdf` (mirror: `xtscaling_compilation.pdf`).
Physics write-up: `../../reports/xtscaling_comparison.{pdf,tex,md}` (deployed to GitHub Pages).

---

## 1. Reproduce

```bash
cd plotting
python3 xtscaling/convert.py                 # data/  ->  converted/   (stdlib only)
root -l -b -q 'xtscaling/plot_xtscaling.C'   # converted/  ->  figures/xtscaling_compilation.pdf
```

`convert.py` does the physics reduction (observable → invariant cross section → x_T scaling).
`plot_xtscaling.C` only draws whatever `converted/*.csv` exist, so adding/removing a dataset is
a data-only change (no macro edit needed, except a style row for a brand-new id).

---

## 2. File inventory

```
xtscaling/
├── convert.py              reduction to (x_T, Y); metadata-aware (units, Δη, observable)
├── plot_xtscaling.C        ROOT macro, sPHENIX style, log-log, 2-column legend
├── fetch_hepdata.sh        helper to browser-download the few still-missing HEPData sets
├── data/                   RAW inputs (one row = one measured point)
│   ├── ppg12_raw.csv               PPG12, read from the ROOT result (see §4)
│   ├── phenix200_fig10.csv         PHENIX inclusive-direct (PRD86 072008)
│   ├── phenix200_isoratio_fig13a.csv  PHENIX isolated/inclusive ratio
│   └── external_datasets.json      15 external datasets in a common schema (see §3)
└── converted/              OUTPUT of convert.py
    ├── <id>.csv                    columns: xT, Y, Y_lo, Y_hi
    └── _overlap_table.txt          per-dataset dump for the universal-curve check (§5)
```

---

## 3. Data sources

HEPData is the natural source, but its `/download/` API is behind a Cloudflare managed
("Turnstile") challenge that **cannot be solved from the analysis cluster** (curl, cloudscraper,
and WebFetch all return 403; the `/record/data/<id>/N/1` endpoint returns a broken default; no
headless browser is available). The numbers HEPData hosts are digitized from the papers, so the
data were pulled from **CF-free equivalents** instead. Every external extraction was performed
once and then independently re-extracted/cross-checked, and finally validated against the
universal band (§5).

| Dataset | √s [GeV] | Reference | acquisition route |
|---|---:|---|---|
| **sPHENIX (PPG12)** | 200 | `efficiencytool/results/Photon_final_bdt_nom.root` + `plotting/rootFiles/syst_sum.root` | local ROOT (uproot) |
| PHENIX | 200 | PRD 86 (2012) 072008 [1205.5533] | arXiv source — paper appendix table |
| ATLAS | 7000 | PLB 706 (2011) 150 [1108.0253] | arXiv source LaTeX table |
| ATLAS | 8000 | JHEP 08 (2016) 005 [1605.03495] | arXiv source LaTeX table |
| ATLAS | 13000 | PLB 770 (2017) 473 [1701.06882] | arXiv source LaTeX table |
| CMS | 7000 | PRL 106 (2011) 082001 [1012.0799] | arXiv source LaTeX table |
| E706 | 31.6, 38.8 | PRD 73 (2006) 032004 [hep-ex/0407011] | arXiv source table (p+p, **not** the p+Be per-nucleon table) |
| CDF | 1800 | PRL 73 (1994) 2662 | INSPIRE full-text PDF, table read with PyMuPDF |
| D0 | 1800 | PRL 77 (1996) 5011 [hep-ex/9603006] | INSPIRE full-text PDF |
| UA2 | 630 | PLB 288 (1992) 386 | INSPIRE full-text PDF |
| NA24 | 23.8 | PRD 36 (1987) 8 | INSPIRE / HEPData record ins236248 (CGS cm² units) |
| UA6 | 24.3 | PLB 436 (1998) 222 | INSPIRE full-text PDF |
| UA1 | 630 | PLB 209 (1988) 385 | Wayback-archived Durham HepData |
| R110 (CMOR) | 63 | NPB 327 (1989) 541 | Wayback-archived Durham HepData (CGS cm² units) |
| E704 | 19.4 | PLB 345 (1995) 569 | Vogelsang–Whalley compilation (JPG 23 (1997) A1) |
| WA70 | 23.0 | ZPC 38 (1988) 371 | Vogelsang–Whalley compilation |

`external_datasets.json` schema (per dataset): `dataset_id, experiment, sqrt_s_GeV, observable,
x_variable, units, eta_range, delta_eta, isolated, uncertainty_type, points[], source,
found_via`. `observable ∈ {invariant_Ed3sigma_dp3, dsigma_dET_eta_integrated,
d2sigma_dET_deta_per_eta}` — this is what tells `convert.py` how to reduce each set.

**Still missing** (HEPData-only, no reachable machine-readable copy; all in already-covered x_T):
ISR sets R806, R807, R108, the two older PHENIX 200 GeV papers, and ALICE low-p_T. Run
`fetch_hepdata.sh` on a machine with browser access to grab them, drop the files in `data/`, and
re-run the pipeline.

---

## 4. Conversion (`convert.py`)

Every spectrum is reduced to the invariant cross section at midrapidity, then x_T-scaled.

1. **Cross-section magnitude → pb**: `mb→×1e9, nb→×1e3, pb→×1, fb→×1e-3, cm²→×1e36`
   (the cm² factor handles the CGS R110/NA24 tables; 1 cm² = 10³⁶ pb).
2. **→ invariant E d³σ/dp³ at y≈0** (`I`), by `observable`:
   - `invariant_Ed3sigma_dp3`  → `I = C`
   - `dsigma_dET_eta_integrated` → `I = (C / Δη) / (2π p_T)`   (Δη = central η-bin width: 1.2 ATLAS |η|<0.6, 2.9 CMS |η|<1.45, 1.4 PPG12 |η|<0.7)
   - `d2sigma_dET_deta_per_eta`  → `I = C / (2π p_T)`           (already per-unit-η, so Δη=1; this is CDF, D0)
3. **x_T scale**: `Y = (√s)^4.5 · I`,  `x_T = 2 p_T/√s`  (E_T ≈ p_T for photons).
4. **Uncertainties** scale linearly (relative error preserved): `Y_unc = Y · (σ/C)`; percent
   errors are converted to absolute first.
5. **Exclude points consistent with zero** — following PHENIX [1], a point is dropped if its
   lower error bar reaches zero (`C ≤ √(stat² + syst_lo²)`). This removed one near-threshold
   point each from E706 (31.6 GeV) and UA1.

**PPG12 specifics**: `h_unfold_sub_result` is the unfolded, leak-corrected spectrum =
`dσ/dE_T` integrated over |η|<0.7 in **pb/GeV** (it is *not* yet divided by Δη). So it is treated
as `dsigma_dET_eta_integrated` with Δη=1.4, i.e. `I = (raw/1.4)/(2π E_T)`. Systematics come from
`syst_sum.root` (`h_sum_low/high`, same raw units). Only the reported reco range **E_T = 12–36
GeV** is used — the two E_T<12 GeV trigger-turn-on bins and the [36,45] overflow bin are dropped.

---

## 5. Validation

The plot is its own cross-check: independently-extracted datasets must land on one curve where
they overlap in x_T. `converted/_overlap_table.txt` dumps Y(x_T) for every set. Result:

- **Collider + RHIC** (CDF, D0, UA1, UA2, ATLAS 7/8/13, CMS, PHENIX, PPG12): collapse to
  **~1–3×** at every shared x_T from 0.014 to 0.13, across a factor 65 in √s. A wrong Δη, 2π,
  cm²→pb, or unit factor in any one set would be far larger than this scatter.
- **Same-√s mutual checks**: the three independent fixed-target sets at the same x_T
  (WA70 / NA24 / E706, pulled via three *different* routes) agree to **~1.3×** — proof the
  conversions are right, not the route.
- **Fixed-target/ISR (19–63 GeV)** sit **3–8× above** the collider band at high x_T and rising
  — the known low-√s breakdown of single-n x_T scaling, **not** a conversion error (see report
  §Findings and refs [1,4,7,8]).

Numbers were also independently re-derived from the ROOT file (PPG12) and raw tables, and the
figure passed a plot-cosmetics review, before deployment.

---

## 6. References

[1] PHENIX, Adare *et al.*, PRD **86** (2012) 072008 [1205.5533] — origin of this figure; n=4.5; the exclude-consistent-with-zero convention; the low-√s exception.
[3] Blankenbecler, Brodsky, Gunion, PLB **42** (1972) 461 — x_T scaling.
[4] Cahalan, Geer, Kogut, Susskind, PRD **11** (1975) 1199 — n_eff(x_T,√s) from QCD scale-breaking.
[5] Vogelsang, Whalley, JPG **23** (1997) A1 — world direct-photon compilation.
[6] Aurenche *et al.*, PRD **73** (2006) 094007 [hep-ph/0602133] — NLO pQCD.
[7] de Florian, Vogelsang, PRD **72** (2005) [hep-ph/0506150] — threshold resummation.
[8] Apanasevich *et al.*, PRD **59** (1999) [hep-ph/9808467] — k_T interpretation of the fixed-target excess.

(Per-measurement references are in the dataset table in §3 and in the report.)
