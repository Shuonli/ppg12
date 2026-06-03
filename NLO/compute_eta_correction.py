#!/usr/bin/env python3
"""
JETPHOX-derived R_eta correction for the PHENIX comparison (PPG12).

R_eta(pT) = (d sigma/(d eta dpT))_{|eta|<0.25} / (d sigma/(d eta dpT))_{|eta|<0.7}
         = (N_|eta|<0.25 / 0.5) / (N_|eta|<0.7 / 1.4)
         = (N_|eta|<0.25 / N_|eta|<0.7) * (1.4 / 0.5)

The denominator is the wider PPG12 acceptance; the numerator is the narrower
PHENIX acceptance. PHENIX points at |eta|<0.25 are divided by R_eta to scale
to |eta|<0.7.

Both numerator and denominator are computed from the SAME chunked CT18NLO
JETPHOX inclusive (no-isolation) run, so PDF/FF/scale dependence cancels in
the ratio. Replaces the previous PYTHIA-derived h_ratio_central_over_full
(efficiencytool/truth_eta_ratio_inclusive.root) for full JETPHOX consistency
with the iso correction (see compute_iso_correction.py).

Input  : rootFiles/jetPHOX_ct18iso200_10_chunked.root (h_truth_eta_pT TH2F,
         eta in [-1, 1] with 0.02-wide bins, pT in [0, 100] with 1 GeV bins).
Output : rootFiles/truth_eta_ratio_jetphox.root (h_ratio_central_over_full).
"""
import os
import sys
import uproot
import numpy as np
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_DEFAULT = os.path.join(HERE, "rootFiles", "jetPHOX_ct18iso200_10_chunked.root")
OUTPUT = os.path.join(HERE, "rootFiles", "truth_eta_ratio_jetphox.root")

input_fn = sys.argv[1] if len(sys.argv) > 1 else INPUT_DEFAULT

f = uproot.open(input_fn)
h2 = f["h_truth_eta_pT"]
values = h2.values()  # shape (n_eta, n_pT)
edges_eta = h2.axes[0].edges()
edges_pt  = h2.axes[1].edges()

# Eta-bin slicing. The TH2F has 0.02-wide bins, so |eta|<0.25 is
# approximated by the bins fully within [-0.24, 0.24] (sub-bin truncation
# is a 4% per-bin geometric effect that mostly cancels in the ratio).
def eta_slice(eta_max):
    lo = np.searchsorted(edges_eta, -eta_max, side='left')
    hi = np.searchsorted(edges_eta,  eta_max, side='right') - 1
    return lo, hi, edges_eta[lo], edges_eta[hi]

lo025, hi025, e025_lo, e025_hi = eta_slice(0.24)  # use 0.24 (closest closed edge inside 0.25)
lo070, hi070, e070_lo, e070_hi = eta_slice(0.70)
dEta_025 = e025_hi - e025_lo
dEta_070 = e070_hi - e070_lo

n_central = values[lo025:hi025, :].sum(axis=0)
n_full    = values[lo070:hi070, :].sum(axis=0)

# Density ratio = (N_narrow / Delta eta_narrow) / (N_wide / Delta eta_wide)
with np.errstate(divide='ignore', invalid='ignore'):
    density_narrow = n_central / dEta_025
    density_wide   = n_full    / dEta_070
    ratio = np.where(density_wide > 0, density_narrow / density_wide, 0.0)

out_f = ROOT.TFile.Open(OUTPUT, "RECREATE")
h_out = ROOT.TH1D(
    "h_ratio_central_over_full",
    "JETPHOX R_eta = density(|eta|<0.25) / density(|eta|<0.7)",
    len(edges_pt) - 1, edges_pt,
)
for i, r in enumerate(ratio):
    h_out.SetBinContent(i + 1, float(r))
h_out.Write()
out_f.Close()

print(f"[compute_eta_correction] input          : {input_fn}")
print(f"[compute_eta_correction] |eta|<0.25 -> eta-bin range [{e025_lo:.3f}, {e025_hi:.3f}], width = {dEta_025:.3f}")
print(f"[compute_eta_correction] |eta|<0.70 -> eta-bin range [{e070_lo:.3f}, {e070_hi:.3f}], width = {dEta_070:.3f}")
print(f"[compute_eta_correction] sample R_eta values:")
for pT_target in (12, 16, 20, 24, 28, 32):
    ib = int(np.searchsorted(edges_pt, pT_target, side='right') - 1)
    if 0 <= ib < len(ratio):
        print(f"  pT = {pT_target:>3d} GeV  -> R_eta = {ratio[ib]:.4f}")
print(f"[compute_eta_correction] wrote {OUTPUT}")
