#!/usr/bin/env python3
"""
Isolation correction for the PHENIX comparison (PPG12).

ratio(pT) = sigma_inclusive(ET_iso=200 GeV) / sigma_isolated(ET_iso=4 GeV, R=0.3)
          = factor to bring an ISOLATED prompt-photon cross section up to the
            INCLUSIVE (no-isolation) definition PHENIX measured.

Numerator : CT18NLO chunked inclusive run produced this session
            (rootFiles/jetPHOX_ct18iso200_10_chunked.root, h_truth_pT)
Denominator: the paper's published isolated CT18NLO curve
            (rootFiles/jetPHOX_ct18_10_chunked.root, h_truth_pT)

Both are d2sigma/(dpT)|_{|y|<0.7}, dir+frag, same PDF/FF/scale/cone, so PDF and
most scale dependence cancel in the ratio. Usage:
    python3 compute_iso_correction.py [incl.root] [iso.root]
"""
import sys
import numpy as np
import uproot

HERE = "/sphenix/user/shuhangli/ppg12/NLO/rootFiles"
incl_fn = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/jetPHOX_ct18iso200_10_chunked.root"
iso_fn  = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/jetPHOX_ct18_10_chunked.root"

def load(fn):
    h = uproot.open(fn)["h_truth_pT"]
    return h.values(), h.errors(), h.axis().edges()

vi, ei, edges_i = load(incl_fn)   # inclusive
vs, es, edges_s = load(iso_fn)    # isolated

if not np.allclose(edges_i, edges_s):
    sys.exit(f"BIN MISMATCH:\n incl={edges_i}\n iso ={edges_s}")
edges = edges_i

# ratio incl/iso and its inverse (iso efficiency), with independent-MC errors
with np.errstate(divide="ignore", invalid="ignore"):
    ratio = vi / vs                                   # >=1 expected
    rel   = np.sqrt((ei/vi)**2 + (es/vs)**2)
    rerr  = ratio * rel
    iso_eff = vs / vi                                 # sigma_iso/sigma_incl, PHENIX measured >0.9
    ieerr   = iso_eff * rel

print(f"\nInclusive : {incl_fn}")
print(f"Isolated  : {iso_fn}\n")
hdr = f"{'pT_lo':>5} {'pT_hi':>5} | {'sig_iso':>10} {'sig_incl':>10} | {'incl/iso':>9} {'+-':>7} | {'iso/incl':>8} {'+-':>6}"
print(hdr); print("-"*len(hdr))
phenix_max = 25.0  # PHENIX 1205.5533 reaches ~25 GeV
for i in range(len(ratio)):
    lo, hi = edges[i], edges[i+1]
    tag = "  <- PHENIX overlap" if hi <= phenix_max + 1e-6 else ""
    print(f"{lo:5.0f} {hi:5.0f} | {vs[i]:10.5g} {vi[i]:10.5g} | "
          f"{ratio[i]:9.4f} {rerr[i]:7.4f} | {iso_eff[i]:8.4f} {ieerr[i]:6.4f}{tag}")

# ---- sanity checks ----
print("\n=== sanity ===")
unphys = [(edges[i], edges[i+1], ratio[i]) for i in range(len(ratio))
          if np.isfinite(ratio[i]) and ratio[i] < 1.0 - 3*rel[i]]
if unphys:
    print(f"  WARNING: {len(unphys)} bin(s) with incl/iso < 1 beyond 3sigma (unphysical):")
    for lo, hi, r in unphys: print(f"    [{lo:.0f},{hi:.0f}] ratio={r:.3f}")
else:
    print("  OK: incl/iso >= 1 within stats in every bin (isolation never increases the rate).")

ov = (edges[:-1] >= 8) & (edges[1:] <= phenix_max + 1e-6)
if ov.any():
    lo_eff = np.nanmin(iso_eff[ov]); hi_eff = np.nanmax(iso_eff[ov])
    print(f"  PHENIX-overlap (8-25 GeV): iso/incl in [{lo_eff:.3f}, {hi_eff:.3f}] "
          f"-> {'consistent' if lo_eff>0.85 else 'TENSION'} with PHENIX measured >0.90")
    print(f"  -> isolation correction (incl/iso) ranges "
          f"{np.nanmin(ratio[ov]):.3f}-{np.nanmax(ratio[ov]):.3f} over the PHENIX overlap")

# ---- save ----
out_txt = f"{HERE}/iso_correction_ct18nlo.txt"
with open(out_txt, "w") as f:
    f.write("# pT_lo pT_hi sig_iso sig_incl incl_over_iso err iso_over_incl err\n")
    for i in range(len(ratio)):
        f.write(f"{edges[i]:.1f} {edges[i+1]:.1f} {vs[i]:.6g} {vi[i]:.6g} "
                f"{ratio[i]:.6f} {rerr[i]:.6f} {iso_eff[i]:.6f} {ieerr[i]:.6f}\n")
print(f"\nwrote {out_txt}")
