#!/usr/bin/env python3
"""
toy_vertex_2d_reweight.py

Toy MC validation of the 2D Gaussian reweight for converting double-interaction
MC (generated with the full beam profile, sigma_MC) to a narrower 1.5 mrad
luminous region (sigma_lum). Demonstrates:

  1. Raw "MC": (z_t1, z_t2) ~ N(0, sigma_MC) x N(0, sigma_MC)
     reco = (z_t1 + z_t2)/2 has sigma_reco = sigma_MC/sqrt(2)

  2. "Target" for 1.5 mrad: both vertices from N(0, sigma_lum)
     reco = (z_t1 + z_t2)/2 has sigma_reco = sigma_lum/sqrt(2)

  3. 2D reweight w(z_t1, z_t2) = [N(0, sigma_lum)/N(0, sigma_MC)]^2:
     After weighting, the toy distributions match the target
     at the percent level.

  4. 1D truth reweight of z_t1 only (keeping z_t2 at full MC width):
     Shows the geometric floor -- reco stays too wide because z_t2
     is not corrected.

Output: reports/figures/vertex_reweight_toy_2d.pdf
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# sPHENIX-style matplotlib rcParams
# ------------------------------------------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
plt.rcParams["mathtext.fontset"] = "stixsans"
plt.rcParams["axes.titlesize"] = 12
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True

np.random.seed(42)

# ------------------------------------------------------------
# Physics parameters (numbers measured from jet12_double at 1.5 mrad)
# ------------------------------------------------------------
SIGMA_MC  = 62.7   # cm, double-MC truth primary vertex width
SIGMA_LUM = 22.3   # cm, 1.5 mrad luminous region (data)
NEVT = 2_000_000   # toy event count

# Draw 2D MC
z_t1_mc = np.random.normal(0.0, SIGMA_MC, NEVT)
z_t2_mc = np.random.normal(0.0, SIGMA_MC, NEVT)
z_r_mc  = 0.5 * (z_t1_mc + z_t2_mc)

# Draw target (correct 1.5 mrad)
z_t1_t  = np.random.normal(0.0, SIGMA_LUM, NEVT)
z_t2_t  = np.random.normal(0.0, SIGMA_LUM, NEVT)
z_r_t   = 0.5 * (z_t1_t + z_t2_t)

# 2D Gaussian reweight: w(z1, z2) = p_target(z1)p_target(z2) / p_MC(z1)p_MC(z2)
# = (sigma_MC^2 / sigma_lum^2) * exp[-0.5*(z1^2+z2^2)*(1/sigma_lum^2 - 1/sigma_MC^2)]
inv_diff = 1.0/SIGMA_LUM**2 - 1.0/SIGMA_MC**2
w_2d = (SIGMA_MC**2 / SIGMA_LUM**2) * np.exp(
    -0.5 * (z_t1_mc**2 + z_t2_mc**2) * inv_diff
)

# 1D truth reweight (keep z_t2 at full MC width): w(z_t1) = p_target(z1)/p_MC(z1)
w_1d_truth = (SIGMA_MC / SIGMA_LUM) * np.exp(
    -0.5 * z_t1_mc**2 * inv_diff
)

# 1D RECO-level reweight (PPG12 nominal production approach):
#   f(z_r) = p_target_reco(z_r) / p_MC_reco(z_r)
# where p_MC_reco is N(0, sigma_MC/sqrt(2)) and p_target_reco is
# N(0, sigma_lum/sqrt(2)). Apply event-by-event using z_r. By
# construction this gives exact reco closure; the question is what it
# does to the truth distributions.
SIGMA_MC_RECO  = SIGMA_MC  / np.sqrt(2.0)
SIGMA_LUM_RECO = SIGMA_LUM / np.sqrt(2.0)
inv_diff_reco = 1.0/SIGMA_LUM_RECO**2 - 1.0/SIGMA_MC_RECO**2
w_1d_reco = (SIGMA_MC_RECO / SIGMA_LUM_RECO) * np.exp(
    -0.5 * z_r_mc**2 * inv_diff_reco
)

# ------------------------------------------------------------
# Build histograms for plotting
# ------------------------------------------------------------
BINS = np.linspace(-200, 200, 81)  # 5 cm bins, same as main analysis
centers = 0.5 * (BINS[1:] + BINS[:-1])

def hist(x, w=None):
    h, _ = np.histogram(x, bins=BINS, weights=w)
    s = h.sum()
    if s > 0:
        h = h / s
    return h

def rms(x, w=None):
    if w is None:
        m = x.mean()
        return float(np.sqrt(((x - m)**2).mean()))
    w = w / w.sum()
    m = (x * w).sum()
    return float(np.sqrt(((x - m)**2 * w).sum()))

hist_z_r_mc   = hist(z_r_mc)
hist_z_r_t    = hist(z_r_t)
hist_z_r_2d   = hist(z_r_mc, w_2d)
hist_z_r_1dt  = hist(z_r_mc, w_1d_truth)
hist_z_r_1dr  = hist(z_r_mc, w_1d_reco)

hist_z_t1_mc  = hist(z_t1_mc)
hist_z_t1_t   = hist(z_t1_t)
hist_z_t1_2d  = hist(z_t1_mc, w_2d)
hist_z_t1_1dt = hist(z_t1_mc, w_1d_truth)
hist_z_t1_1dr = hist(z_t1_mc, w_1d_reco)

hist_z_t2_mc  = hist(z_t2_mc)
hist_z_t2_t   = hist(z_t2_t)
hist_z_t2_2d  = hist(z_t2_mc, w_2d)
hist_z_t2_1dt = hist(z_t2_mc, w_1d_truth)
hist_z_t2_1dr = hist(z_t2_mc, w_1d_reco)

print("Width summary (cm RMS):")
print(f"  MC raw:                 z_t1={rms(z_t1_mc):.2f}  z_t2={rms(z_t2_mc):.2f}  z_r={rms(z_r_mc):.2f}")
print(f"  Target (direct):        z_t1={rms(z_t1_t):.2f}   z_t2={rms(z_t2_t):.2f}   z_r={rms(z_r_t):.2f}")
print(f"  MC x f(z_r) [reco]:     z_t1={rms(z_t1_mc, w_1d_reco):.2f}  z_t2={rms(z_t2_mc, w_1d_reco):.2f}  z_r={rms(z_r_mc, w_1d_reco):.2f}")
print(f"  MC x w(z_t1) [truth 1D]: z_t1={rms(z_t1_mc, w_1d_truth):.2f}  z_t2={rms(z_t2_mc, w_1d_truth):.2f}  z_r={rms(z_r_mc, w_1d_truth):.2f}")
print(f"  MC x w(z_t1,z_t2) [2D]: z_t1={rms(z_t1_mc, w_2d):.2f}  z_t2={rms(z_t2_mc, w_2d):.2f}  z_r={rms(z_r_mc, w_2d):.2f}")

# TV distance on |z|<60 (matching the main analysis vertex cut)
mask60 = np.abs(centers) < 60
def tv(h1, h2):
    return 0.5 * np.abs(h1[mask60] - h2[mask60]).sum()

def tv_all(h_mc_z_t1, h_mc_z_t2, h_mc_z_r, tag):
    t1 = tv(h_mc_z_t1, hist_z_t1_t)
    t2 = tv(h_mc_z_t2, hist_z_t2_t)
    tr = tv(h_mc_z_r,  hist_z_r_t)
    print(f"  {tag:<28}  z_t1={t1:.4f}   z_t2={t2:.4f}   z_r={tr:.4f}")

print("\nTV distance on |z|<60 between reweighted MC and target, per axis:")
tv_all(hist_z_t1_mc,  hist_z_t2_mc,  hist_z_r_mc,  "raw MC")
tv_all(hist_z_t1_1dr, hist_z_t2_1dr, hist_z_r_1dr, "MC x f(z_r)     [reco]")
tv_all(hist_z_t1_1dt, hist_z_t2_1dt, hist_z_r_1dt, "MC x w(z_t1)   [truth 1D]")
tv_all(hist_z_t1_2d,  hist_z_t2_2d,  hist_z_r_2d,  "MC x w(z_t1,z_t2) [2D]")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

def panel(ax, title, hD_mc, hD_t, hD_1dr, hD_1dt, hD_2d, xlabel):
    ax.step(centers, hD_mc, where="mid",
            label="raw MC (full profile)", color="0.5", linewidth=2, linestyle="-")
    ax.step(centers, hD_t,  where="mid",
            label="target (correct 1.5 mrad)", color="k", linewidth=2.5, linestyle="-")
    ax.step(centers, hD_1dr, where="mid",
            label=r"MC $\times\,f(z_{r})$  (nominal production, reco-level)",
            color="C2", linewidth=2, linestyle="-.")
    ax.step(centers, hD_1dt, where="mid",
            label=r"MC $\times\,w(z_{t,1})$  (truth 1D only)",
            color="C3", linewidth=2, linestyle=":")
    ax.step(centers, hD_2d, where="mid",
            label=r"MC $\times\,w(z_{t,1},z_{t,2})$  (2D truth)",
            color="C0", linewidth=2, linestyle="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("normalized counts")
    ax.set_xlim(-120, 120)
    # Clear headroom above the tallest curve drawn in this panel
    max_observed_peak = max(hD_mc.max(), hD_t.max(),
                            hD_1dr.max(), hD_1dt.max(), hD_2d.max())
    ax.set_ylim(0, 1.50 * max_observed_peak)
    ax.set_title(title, fontsize=12)
    ax.legend(loc="upper right", fontsize=9, frameon=False,
              handlelength=1.8, labelspacing=0.35)
    ax.axvline(-60, ls=":", color="0.7", lw=1)
    ax.axvline(+60, ls=":", color="0.7", lw=1)
    # sPHENIX branding (upper-left corner)
    ax.text(0.04, 0.96, r"$\bf{\it{sPHENIX}}$ Internal",
            transform=ax.transAxes, fontsize=11, va="top", ha="left")
    ax.text(0.04, 0.89, r"$p+p$  $\sqrt{s}$ = 200 GeV",
            transform=ax.transAxes, fontsize=9, va="top", ha="left")
    return max_observed_peak

peak_t1 = panel(axes[0,0], r"Primary truth vertex  $z_{t,1}$",
      hist_z_t1_mc, hist_z_t1_t, hist_z_t1_1dr, hist_z_t1_1dt, hist_z_t1_2d,
      r"$z_{t,1}$ [cm]")
peak_t2 = panel(axes[0,1], r"Secondary truth vertex  $z_{t,2}$",
      hist_z_t2_mc, hist_z_t2_t, hist_z_t2_1dr, hist_z_t2_1dt, hist_z_t2_2d,
      r"$z_{t,2}$ [cm]")
peak_zr = panel(axes[1,0], r"Reco vertex  $z_{r}=(z_{t,1}+z_{t,2})/2$",
      hist_z_r_mc, hist_z_r_t, hist_z_r_1dr, hist_z_r_1dt, hist_z_r_2d,
      r"$z_{r}$ [cm]")
print("\nPer-panel peak heights and applied y-limits (y_max = 1.50 * peak):")
print(f"  z_t,1 panel: peak={peak_t1:.4f}  ylim=(0, {1.50*peak_t1:.4f})")
print(f"  z_t,2 panel: peak={peak_t2:.4f}  ylim=(0, {1.50*peak_t2:.4f})")
print(f"  z_r   panel: peak={peak_zr:.4f}  ylim=(0, {1.50*peak_zr:.4f})")

# Panel 4: RMS widths + TV distances as a compact table
ax = axes[1,1]
ax.axis("off")

# Header lines use the regular sans-serif family; the numerical columns
# below are rendered with a monospace font so they line up cleanly.
header_lines = []
header_lines.append((r"$\bf{Toy\ MC\ width\ summary}$  (RMS, cm)", "sans"))
header_lines.append(("", "sans"))
header_lines.append((r"$\sigma_{\mathrm{MC}} = 62.7$ cm  (double-MC beam profile)", "sans"))
header_lines.append((r"$\sigma_{\mathrm{lum}} = 22.3$ cm  (1.5 mrad data)", "sans"))
header_lines.append(("", "sans"))
header_lines.append((r"$\bf{RMS\ widths}$   ($z_{t,1}$,  $z_{t,2}$,  $z_r$):", "sans"))

# Monospace rows: fixed-width label column + fixed-width number columns.
def _row(label, a, b, c):
    return f"  {label:<22}{a:>6.1f}  {b:>6.1f}  {c:>6.1f}"

mono_lines = []
mono_lines.append(_row("raw MC",           rms(z_t1_mc),           rms(z_t2_mc),           rms(z_r_mc)))
mono_lines.append(_row("target",           rms(z_t1_t),            rms(z_t2_t),            rms(z_r_t)))
mono_lines.append(_row("f(z_r) reco",      rms(z_t1_mc, w_1d_reco), rms(z_t2_mc, w_1d_reco), rms(z_r_mc, w_1d_reco)))
mono_lines.append(_row("w(z_t1) truth 1D", rms(z_t1_mc, w_1d_truth), rms(z_t2_mc, w_1d_truth), rms(z_r_mc, w_1d_truth)))
mono_lines.append(_row("w(z_t1,z_t2) 2D",  rms(z_t1_mc, w_2d),      rms(z_t2_mc, w_2d),      rms(z_r_mc, w_2d)))

tv_header = (r"$\bf{TV\ distance\ |z|<60}$  ($z_r$ axis):", "sans")

def _tv_row(label, val, note):
    return f"  {label:<22}{val:>7.4f}   {note}"

tv_lines = []
tv_lines.append(_tv_row("raw MC",           tv(hist_z_r_mc,  hist_z_r_t), ""))
tv_lines.append(_tv_row("f(z_r) reco",      tv(hist_z_r_1dr, hist_z_r_t), "(closes reco)"))
tv_lines.append(_tv_row("w(z_t1) truth 1D", tv(hist_z_r_1dt, hist_z_r_t), "(fails)"))
tv_lines.append(_tv_row("w(z_t1,z_t2) 2D",  tv(hist_z_r_2d,  hist_z_r_t), "(closes everything)"))

# Line step: ~0.048 fits 20+ lines comfortably; start near the top.
y = 0.97
dy = 0.048
for text, kind in header_lines:
    ax.text(0.02, y, text, transform=ax.transAxes, fontsize=10, va="top")
    y -= dy
for text in mono_lines:
    ax.text(0.02, y, text, transform=ax.transAxes,
            fontsize=10, family="monospace", va="top")
    y -= dy
y -= 0.02
ax.text(0.02, y, tv_header[0], transform=ax.transAxes, fontsize=10, va="top")
y -= dy
for text in tv_lines:
    ax.text(0.02, y, text, transform=ax.transAxes,
            fontsize=10, family="monospace", va="top")
    y -= dy

fig.suptitle(
    r"Toy MC: 2D vs 1D truth-vertex reweight of double-interaction MC",
    fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])

import os
os.makedirs("reports/figures", exist_ok=True)
out = "reports/figures/vertex_reweight_toy_2d.pdf"
fig.savefig(out)
print(f"\nWrote {out}")
