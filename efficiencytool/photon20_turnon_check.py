#!/usr/bin/env python3
# Check if the photon10/photon20 truth-pT boundary can be lowered from 30 GeV.
# Overlays the unfiltered d sigma/dpT spectra from photon5, photon10, photon20
# to locate (a) photon20's pT-hat turn-on and (b) the pT range where photon20
# and photon10 agree in the weighted differential cross section.

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

RESULTS = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
FIG = "/sphenix/user/shuhangli/ppg12/plotting/figures/truth_spectrum"

fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10, 9),
                              gridspec_kw={"height_ratios": [3, 1]},
                              sharex=True)

samples = [("photon5",  "tab:red",   "photon5  (pT-hat > 5)"),
           ("photon10", "tab:blue",  "photon10 (pT-hat > 10)"),
           ("photon20", "tab:green", "photon20 (pT-hat > 20)")]

ref_vals = None
ref_errs = None
ref_centers = None

plot_data = {}
for (s, color, label) in samples:
    f = uproot.open(f"{RESULTS}/truth_spectrum_{s}.root")
    h = f["h_max_truth_photon_pt"]
    vals = h.values()
    errs = h.errors()
    edges = h.axes[0].edges()
    centers = 0.5 * (edges[:-1] + edges[1:])
    per_w = float(str(f["per_entry_weight_pb_per_GeV"].member("fTitle")))
    n_entries = np.round(vals / per_w).astype(int) if per_w > 0 else vals * 0
    stat_err = per_w * np.sqrt(np.maximum(n_entries, 0))
    plot_data[s] = {"centers": centers, "vals": vals, "err": stat_err,
                    "n_entries": n_entries, "label": label, "color": color}

# Top: overlay of d sigma/dpT
for s, d in plot_data.items():
    m = (d["centers"] >= 14) & (d["centers"] <= 40)
    ax.errorbar(d["centers"][m], d["vals"][m], yerr=d["err"][m],
                fmt="o", markersize=3, color=d["color"],
                label=d["label"], alpha=0.85)

ax.set_yscale("log")
ax.set_ylabel(r"d$\sigma$/d$p_T$ [pb/GeV]")
ax.set_title("photon5/10/20 unfiltered truth-photon spectrum — overlap region 14-40 GeV")

# Mark pT-hat thresholds and boundaries
ax.axvline(20, color="gray", linestyle=":", alpha=0.5)
ax.text(20.1, 5e2, "photon20 pT-hat", rotation=90, fontsize=8, va="top", color="tab:green")
ax.axvline(30, color="red", linestyle="--", alpha=0.7)
ax.text(30.1, 5e2, "current p20\nlower bound", rotation=90, fontsize=8, va="top", color="red")

# Candidate boundaries to evaluate
for p, lbl in [(22, "candidate: 22"), (24, "candidate: 24"), (25, "candidate: 25")]:
    ax.axvline(p, color="purple", linestyle=":", alpha=0.4)

ax.legend(loc="upper right", frameon=False, fontsize=10)
ax.grid(True, which="both", alpha=0.3)

# Bottom: ratio photon20/photon10 (or photon5/photon10 as reference) to see where they agree
m20 = plot_data["photon20"]["vals"] > 0
m10 = plot_data["photon10"]["vals"] > 0
m = m10 & m20
r20_10 = plot_data["photon20"]["vals"][m] / plot_data["photon10"]["vals"][m]
# Errors
v20 = plot_data["photon20"]["vals"][m]
v10 = plot_data["photon10"]["vals"][m]
e20 = plot_data["photon20"]["err"][m]
e10 = plot_data["photon10"]["err"][m]
rel_err = np.sqrt((e20 / v20) ** 2 + (e10 / v10) ** 2)
r_err = r20_10 * rel_err
centers_m = plot_data["photon20"]["centers"][m]
mm = (centers_m >= 14) & (centers_m <= 40)

ax2.errorbar(centers_m[mm], r20_10[mm], yerr=r_err[mm],
             fmt="o", markersize=3, color="tab:green",
             label=r"photon20 / photon10")
ax2.axhline(1.0, color="black", linestyle="--", alpha=0.6)
ax2.axhspan(0.95, 1.05, alpha=0.15, color="green", label=r"$\pm 5\%$ agreement")
ax2.set_ylim(0.0, 2.0)
ax2.set_xlim(14, 40)
ax2.set_xlabel("max truth photon pT [GeV]")
ax2.set_ylabel("ratio")
ax2.axvline(20, color="gray", linestyle=":", alpha=0.5)
ax2.axvline(30, color="red", linestyle="--", alpha=0.7)
ax2.axvline(22, color="purple", linestyle=":", alpha=0.4)
ax2.axvline(24, color="purple", linestyle=":", alpha=0.4)
ax2.axvline(25, color="purple", linestyle=":", alpha=0.4)
ax2.legend(loc="upper right", frameon=False, fontsize=9)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
outpath = f"{FIG}/photon20_turnon_check.pdf"
fig.savefig(outpath)
plt.close(fig)
print(f"wrote {outpath}")

# Print stats comparison
print()
print("Statistics comparison (per 0.6 GeV bin):")
print(f"{'pT (GeV)':>10} | {'p10 N':>10} | {'p10 rel σ':>10} | {'p20 N':>10} | {'p20 rel σ':>10} | {'p20/p10':>8} | {'σ gain':>8}")
for pT_test in [20, 22, 24, 25, 26, 28, 30]:
    idx = np.argmin(np.abs(plot_data["photon10"]["centers"] - pT_test))
    n10 = plot_data["photon10"]["n_entries"][idx]
    n20 = plot_data["photon20"]["n_entries"][idx]
    v10 = plot_data["photon10"]["vals"][idx]
    v20 = plot_data["photon20"]["vals"][idx]
    e10 = plot_data["photon10"]["err"][idx]
    e20 = plot_data["photon20"]["err"][idx]
    r10 = e10 / v10 if v10 > 0 else 0
    r20 = e20 / v20 if v20 > 0 else 0
    ratio = v20 / v10 if v10 > 0 else 0
    gain = r10 / r20 if r20 > 0 else 0
    print(f"{pT_test:>10} | {n10:>10} | {r10:>9.2%} | {n20:>10} | {r20:>9.2%} | {ratio:>8.3f} | {gain:>7.1f}x")
