#!/usr/bin/env python3
"""JETPHOX NLO statistical-precision study for PPG12.

Computes:
  1. Per-bin JETPHOX stat_rel for mu=0.5/1.0/2.0 x ET scales.
  2. Per-bin measured total uncertainty from Photon_final_bdt_nom.root +
     syst_sum.root (stat from unfolded histo, syst from h_sum_rel_{high,low}).
  3. Ratio JETPHOX stat_rel / measured total_rel per bin, flagging bins where
     JETPHOX stat is >= T/3 (i.e. theory stat not much smaller than the total
     measurement uncertainty).
  4. Scale-envelope relative uncertainty vs JETPHOX stat_rel (honesty check
     for the scale band).
  5. Event-count scaling factors needed in each bin to reach stat_rel = T/3.

Outputs three CSV files plus a markdown summary.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import uproot

BASE = Path("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12")
NLO_DIR = BASE / "NLO" / "rootFiles"
EFF_RESULTS = BASE / "efficiencytool" / "results"
SYST_FILE = BASE / "plotting" / "rootFiles" / "syst_sum.root"
OUT_DIR = BASE / "reports" / "work_jetphox_stats"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Fiducial |eta| < 0.7 -> delta_eta = 1.4; applied in plot_final.C.
DETA = 1.4

# Target JETPHOX stat_rel relative to measured total uncertainty.
# If JETPHOX stat_rel <= TARGET_FRAC * total_rel then theory stat is negligible.
TARGET_FRAC = 1.0 / 3.0


def read_jetphox(scale: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (edges, values, errors) of h_truth_pT (raw, unscaled by deta)."""
    f = uproot.open(NLO_DIR / f"jetPHOX_{scale}.root")
    h = f["h_truth_pT"]
    edges = h.axis().edges()
    return edges, h.values(), h.errors()


def read_measured() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (edges, values, stat errors) of h_unfold_sub_result (raw)."""
    f = uproot.open(EFF_RESULTS / "Photon_final_bdt_nom.root")
    h = f["h_unfold_sub_result"]
    edges = h.axis().edges()
    return edges, h.values(), h.errors()


def read_syst_rel() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (edges, rel_low, rel_high) of the total systematic band."""
    f = uproot.open(SYST_FILE)
    h_low = f["h_sum_rel_low"]
    h_high = f["h_sum_rel_high"]
    edges = h_low.axis().edges()
    return edges, h_low.values(), h_high.values()


def bin_centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def main() -> None:
    edges_j, v05, e05 = read_jetphox("05")
    _, v10, e10 = read_jetphox("10")
    _, v20, e20 = read_jetphox("20")

    edges_m, meas_v, meas_stat = read_measured()
    edges_s, syst_low, syst_high = read_syst_rel()

    # All three sets share the same 14-bin truth binning [7, 8, 10, ..., 45].
    assert np.allclose(edges_j, edges_m)
    assert np.allclose(edges_j, edges_s)

    edges = edges_j
    lows = edges[:-1]
    highs = edges[1:]
    widths = highs - lows
    centers = bin_centers(edges)

    # Apply the same 1/deta normalization the plot applies -- does not affect
    # relative uncertainties, but makes absolute numbers match plot-ready
    # d^2 sigma / d eta / d pT [pb/GeV].
    j05 = v05 / DETA
    j10 = v10 / DETA
    j20 = v20 / DETA
    j05_err = e05 / DETA
    j10_err = e10 / DETA
    j20_err = e20 / DETA

    meas = meas_v / DETA
    meas_err = meas_stat / DETA

    # JETPHOX "number of effective events per bin" estimator.
    # If sum of weights (value) = V and sum of weights squared = E^2,
    # then N_eff = V^2 / E^2. Use RAW (pre-/deta) values because /deta is a
    # constant global factor that cancels in the ratio.
    def neff(v, e):
        out = np.zeros_like(v)
        mask = (v > 0) & (e > 0)
        out[mask] = (v[mask] / e[mask]) ** 2
        return out

    neff05 = neff(v05, e05)
    neff10 = neff(v10, e10)
    neff20 = neff(v20, e20)

    # --- Table 1: per-bin JETPHOX stats (all three scales) ------------------
    t1_path = OUT_DIR / "table1_jetphox_per_bin_stats.csv"
    with t1_path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "scale",
                "pT_low",
                "pT_high",
                "pT_center",
                "dsigma_dpT_pb_per_GeV",
                "stat_abs_pb_per_GeV",
                "stat_rel_pct",
                "N_eff_in_bin",
            ]
        )
        for scale, val, err, neff_arr in [
            ("05", j05, j05_err, neff05),
            ("10", j10, j10_err, neff10),
            ("20", j20, j20_err, neff20),
        ]:
            for i in range(len(val)):
                rel = (err[i] / val[i] * 100.0) if val[i] > 0 else float("nan")
                w.writerow(
                    [
                        scale,
                        f"{lows[i]:.1f}",
                        f"{highs[i]:.1f}",
                        f"{centers[i]:.2f}",
                        f"{val[i]:.6g}",
                        f"{err[i]:.6g}",
                        f"{rel:.3f}",
                        f"{neff_arr[i]:.2f}",
                    ]
                )

    # --- Table 2: data vs theory precision ----------------------------------
    # Measured total_rel = sqrt(stat_rel^2 + syst_rel^2), with syst_rel =
    # max(|rel_low|, |rel_high|) for a symmetric conservative proxy.
    j10_rel = np.where(j10 > 0, j10_err / j10, np.nan)
    meas_rel_stat = np.where(meas > 0, meas_err / meas, np.nan)

    syst_rel_max = np.maximum(np.abs(syst_low), np.abs(syst_high))
    meas_rel_total = np.sqrt(meas_rel_stat**2 + syst_rel_max**2)

    ratio = j10_rel / meas_rel_total

    # Scale factor on JETPHOX event count required to reach stat_rel=TARGET_FRAC*total.
    target = TARGET_FRAC * meas_rel_total
    scale_factor = np.where(target > 0, (j10_rel / target) ** 2, np.nan)

    t2_path = OUT_DIR / "table2_data_vs_theory_precision.csv"
    with t2_path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "pT_low",
                "pT_high",
                "pT_center",
                "dsigma_meas_pb_per_GeV",
                "meas_stat_rel_pct",
                "meas_syst_rel_high_pct",
                "meas_syst_rel_low_pct",
                "meas_total_rel_pct",
                "jetphox_stat_rel_pct",
                "ratio_jstat_over_total",
                "verdict",
                "N_scale_to_reach_T_over_3",
            ]
        )
        for i in range(len(meas)):
            verdict = (
                "OK"
                if (j10_rel[i] <= TARGET_FRAC * meas_rel_total[i])
                else ("NEEDS MORE" if j10_rel[i] > 0.5 * meas_rel_total[i] else "MARGINAL")
            )
            w.writerow(
                [
                    f"{lows[i]:.1f}",
                    f"{highs[i]:.1f}",
                    f"{centers[i]:.2f}",
                    f"{meas[i]:.6g}",
                    f"{meas_rel_stat[i]*100:.3f}",
                    f"{syst_high[i]*100:.3f}",
                    f"{syst_low[i]*100:.3f}",
                    f"{meas_rel_total[i]*100:.3f}",
                    f"{j10_rel[i]*100:.3f}",
                    f"{ratio[i]:.3f}",
                    verdict,
                    f"{scale_factor[i]:.3f}",
                ]
            )

    # --- Table 3: scale envelope vs stat ------------------------------------
    env_abs = 0.5 * (np.maximum(j05, j20) - np.minimum(j05, j20))
    env_rel = np.where(j10 > 0, env_abs / j10, np.nan)
    env_over_stat = env_rel / j10_rel

    t3_path = OUT_DIR / "table3_envelope_vs_stat.csv"
    with t3_path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "pT_low",
                "pT_high",
                "pT_center",
                "scale_envelope_rel_pct",
                "jetphox_stat_rel_pct",
                "envelope_over_stat_ratio",
                "theory_band_honest",
            ]
        )
        for i in range(len(j10)):
            honest = env_rel[i] >= 3.0 * j10_rel[i]
            w.writerow(
                [
                    f"{lows[i]:.1f}",
                    f"{highs[i]:.1f}",
                    f"{centers[i]:.2f}",
                    f"{env_rel[i]*100:.3f}",
                    f"{j10_rel[i]*100:.3f}",
                    f"{env_over_stat[i]:.3f}",
                    "YES" if honest else "NO-stat-contaminates-band",
                ]
            )

    print("Wrote", t1_path)
    print("Wrote", t2_path)
    print("Wrote", t3_path)

    # --- Stdout summary ----------------------------------------------------
    print()
    print("=== Per-bin JETPHOX stat_rel summary (nominal mu=1.0) ===")
    print(
        f"{'bin':>14}  {'d sigma/dpT':>12}  {'stat_rel%':>10}  {'N_eff':>10}"
    )
    for i in range(len(j10)):
        print(
            f"  [{lows[i]:>4.1f}, {highs[i]:>4.1f}]  {j10[i]:12.4g}  "
            f"{j10_rel[i]*100:10.3f}  {neff10[i]:10.1f}"
        )

    print()
    print("=== Data vs theory precision ===")
    print(
        f"{'bin':>14}  {'meas':>10}  {'stat%':>7}  {'syst%':>7}  "
        f"{'tot%':>7}  {'jstat%':>7}  {'ratio':>6}  verdict"
    )
    for i in range(len(meas)):
        verdict = (
            "OK"
            if (j10_rel[i] <= TARGET_FRAC * meas_rel_total[i])
            else ("NEEDS" if j10_rel[i] > 0.5 * meas_rel_total[i] else "MARGINAL")
        )
        print(
            f"  [{lows[i]:>4.1f}, {highs[i]:>4.1f}]  {meas[i]:10.4g}  "
            f"{meas_rel_stat[i]*100:7.2f}  {syst_rel_max[i]*100:7.2f}  "
            f"{meas_rel_total[i]*100:7.2f}  {j10_rel[i]*100:7.2f}  "
            f"{ratio[i]:6.2f}  {verdict}"
        )

    print()
    print("=== Scale envelope vs stat ===")
    print(
        f"{'bin':>14}  {'env%':>7}  {'jstat%':>7}  {'env/stat':>9}  honest?"
    )
    for i in range(len(j10)):
        honest = env_rel[i] >= 3.0 * j10_rel[i]
        print(
            f"  [{lows[i]:>4.1f}, {highs[i]:>4.1f}]  "
            f"{env_rel[i]*100:7.2f}  {j10_rel[i]*100:7.2f}  "
            f"{env_over_stat[i]:9.2f}  {'YES' if honest else 'NO'}"
        )


if __name__ == "__main__":
    main()
