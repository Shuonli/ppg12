#!/usr/bin/env python3
"""Make converted-vs-unconverted photon response plots from response_converted.root.

Produces:
  figures/response_shape_grid.pdf       — R_E per-pT panel, categories overlaid
  figures/response_shape_ET_grid.pdf    — R_ET per-pT panel
  figures/response_mean_vs_pt.pdf       — <R_E> and <R_ET> vs pT
  figures/response_sigma_vs_pt.pdf      — Gaussian core sigma(R_E) vs pT
  figures/response_tailfrac_vs_pt.pdf   — frac R<0.8 and R<0.5 vs pT
  figures/response_angular_deltas.pdf   — dEta, dPhi, dR per category
  figures/response_match_eff_vs_pt.pdf  — matching efficiency vs pT
  figures/response_eta_dep.pdf          — mean R_E and sigma vs eta
  rootFiles/summary.pkl                 — numeric summary arrays
"""
from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
ROOT_PATH = HERE / "rootFiles" / "response_converted.root"
SAMPLE_PATH = HERE / "rootFiles" / "samples.pkl"
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PT_BINS = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
ETA_BINS = np.array([-0.7, -0.35, 0.0, 0.35, 0.7], dtype=float)
R_MIN, R_MAX, R_NBIN = 0.0, 1.3, 130

CATS = [
    ("unconv", "Unconverted", "k", "o", "-"),
    ("conv", "Converted", "r", "s", "--"),
    ("bad", "Bad-secondary", "b", "^", ":"),
]
CAT_NAMES = [c[0] for c in CATS]


def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def crystal_ball(x, N, mu, sigma, alpha, n):
    """Single-sided Crystal Ball function with left (low-R) power-law tail.

    Standard calorimeter parametrisation: Gaussian core for (x-mu)/sigma > -alpha,
    power-law tail for (x-mu)/sigma <= -alpha.
    """
    sigma = abs(sigma)
    alpha = abs(alpha)
    t = (x - mu) / sigma
    A_coef = (n / alpha) ** n * np.exp(-0.5 * alpha * alpha)
    B_coef = n / alpha - alpha
    # Guard: B - t > 0 is required; for t >> -alpha the tail branch is not used,
    # so we clip before the power.
    tail_base = np.clip(B_coef - t, 1e-8, None)
    result = np.where(t > -alpha,
                      np.exp(-0.5 * t * t),
                      A_coef * np.power(tail_base, -n))
    return N * result


def read_hist(f, name):
    h = f[name]
    counts, edges = h.to_numpy()
    centers = 0.5 * (edges[:-1] + edges[1:])
    return counts, edges, centers


def read_2d(f, name):
    h = f[name]
    c, xe, ye = h.to_numpy()
    return c, xe, ye


def _find_peak(counts, centers):
    """Return (peak_idx, peak_x, peak_y) via 3-bin smoothed argmax."""
    smooth = np.convolve(counts, np.ones(3) / 3.0, mode="same")
    idx = int(np.argmax(smooth))
    return idx, float(centers[idx]), float(counts[idx])


def gauss_fit(counts, centers, min_counts=50, half_width=0.10):
    """Single-Gaussian fit in a symmetric window of half-width 0.10 around the peak.

    Fits the core of the R distribution without the low-R tail — the result is
    the 'core' Gaussian mu and sigma comparable to a standard calorimeter
    resolution fit. Returns None on failure or clearly pathological result.
    """
    if counts.sum() < min_counts:
        return None
    peak_idx, peak_x, peak_y = _find_peak(counts, centers)
    lo = peak_x - half_width
    hi = peak_x + half_width
    msk = (centers >= lo) & (centers <= hi)
    if msk.sum() < 5 or counts[msk].sum() < min_counts:
        return None
    x = centers[msk]
    y = counts[msk]
    y_err = np.sqrt(np.maximum(y, 1.0))
    try:
        popt, pcov = curve_fit(
            gauss, x, y,
            p0=[peak_y, peak_x, 0.07],
            sigma=y_err, absolute_sigma=False,
            maxfev=5000,
        )
    except Exception:
        return None
    perr = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else [np.nan] * 3
    if abs(popt[1] - peak_x) > 0.15 or abs(popt[2]) > 0.3 or abs(popt[2]) < 1e-4:
        return None
    return {"A": popt[0], "mu": popt[1], "sigma": abs(popt[2]),
            "A_err": perr[0], "mu_err": perr[1], "sigma_err": perr[2],
            "lo": lo, "hi": hi}


def cb_fit(counts, centers, min_counts=500, fit_lo=0.3, fit_hi=1.2):
    """Fit a Crystal Ball function (left-side power-law tail) on the response.

    Captures both the Gaussian core and the asymmetric low-R tail caused by
    energy leakage / conversions. Returns None on fit failure.
    Parameters: N (norm), mu, sigma (core), alpha (tail transition in units of
    sigma), n (tail power-law index).
    """
    if counts.sum() < min_counts:
        return None
    peak_idx, peak_x, peak_y = _find_peak(counts, centers)
    msk = (centers >= fit_lo) & (centers <= fit_hi)
    if msk.sum() < 20:
        return None
    x = centers[msk]
    y = counts[msk]
    y_err = np.sqrt(np.maximum(y, 1.0))
    p0 = [peak_y, max(min(peak_x, 1.05), 0.5), 0.07, 1.5, 5.0]
    bounds = ([0.0, 0.4, 0.005, 0.1, 1.01],
              [np.inf, 1.1, 0.3, 10.0, 200.0])
    try:
        popt, pcov = curve_fit(
            crystal_ball, x, y,
            p0=p0, bounds=bounds,
            sigma=y_err, absolute_sigma=False,
            maxfev=20000,
        )
    except Exception:
        return None
    perr = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else [np.nan] * 5
    if abs(popt[1] - peak_x) > 0.2 or abs(popt[2]) < 1e-4:
        return None
    return {"N": popt[0], "mu": popt[1], "sigma": abs(popt[2]),
            "alpha": abs(popt[3]), "n": popt[4],
            "N_err": perr[0], "mu_err": perr[1], "sigma_err": perr[2],
            "alpha_err": perr[3], "n_err": perr[4],
            "lo": fit_lo, "hi": fit_hi}


def weighted_moments(counts, centers):
    tot = counts.sum()
    if tot < 1:
        return np.nan, np.nan
    mu = np.sum(centers * counts) / tot
    var = np.sum((centers - mu) ** 2 * counts) / tot
    return mu, np.sqrt(max(var, 0.0))


def weighted_median_percentiles(counts, centers, qs=(0.16, 0.50, 0.84)):
    tot = counts.sum()
    if tot < 1:
        return [np.nan] * len(qs)
    cdf = np.cumsum(counts) / tot
    out = []
    for q in qs:
        idx = np.searchsorted(cdf, q, side="left")
        idx = min(idx, len(centers) - 1)
        out.append(centers[idx])
    return out


def tail_fraction(counts, centers, threshold):
    tot = counts.sum()
    if tot < 1:
        return np.nan, np.nan
    msk = centers < threshold
    n = counts[msk].sum()
    frac = n / tot
    err = np.sqrt(frac * (1 - frac) / tot) if tot > 0 else np.nan
    return frac, err


def sphenix_label(ax, x=0.035, y=0.96, extra=None):
    ax.text(x, y, r"$\bf{sPHENIX}$ Simulation", transform=ax.transAxes,
            ha="left", va="top", fontsize=10)
    ax.text(x, y - 0.065, r"photon10, $|\eta|<0.7$",
            transform=ax.transAxes, ha="left", va="top", fontsize=8, color="gray")
    if extra:
        ax.text(x, y - 0.125, extra, transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color="gray")


def make_summary(f):
    """Compute per-(pT, cat) summary: moments, percentiles, Gaussian + CB fits."""
    out = {c: {
        "mean": [], "rms": [],
        "gauss_mu": [], "gauss_sigma": [],
        "cb_mu": [], "cb_sigma": [], "cb_alpha": [], "cb_n": [],
        "median": [], "p16": [], "p84": [],
        "tail_0p8": [], "tail_0p5": [],
        "n": [],
    } for c in CAT_NAMES}
    R_edges = np.linspace(R_MIN, R_MAX, R_NBIN + 1)
    centers = 0.5 * (R_edges[:-1] + R_edges[1:])
    for ip in range(len(PT_BINS) - 1):
        for cat in CAT_NAMES:
            counts, _, _ = read_hist(f, f"R_E_{cat}_pt{ip}")
            mu, rms = weighted_moments(counts, centers)
            p16, p50, p84 = weighted_median_percentiles(counts, centers)
            t8, _ = tail_fraction(counts, centers, 0.8)
            t5, _ = tail_fraction(counts, centers, 0.5)
            g = gauss_fit(counts, centers)
            cb = cb_fit(counts, centers)
            out[cat]["mean"].append(mu)
            out[cat]["rms"].append(rms)
            out[cat]["gauss_mu"].append(g["mu"] if g else np.nan)
            out[cat]["gauss_sigma"].append(g["sigma"] if g else np.nan)
            out[cat]["cb_mu"].append(cb["mu"] if cb else np.nan)
            out[cat]["cb_sigma"].append(cb["sigma"] if cb else np.nan)
            out[cat]["cb_alpha"].append(cb["alpha"] if cb else np.nan)
            out[cat]["cb_n"].append(cb["n"] if cb else np.nan)
            out[cat]["median"].append(p50)
            out[cat]["p16"].append(p16)
            out[cat]["p84"].append(p84)
            out[cat]["tail_0p8"].append(t8)
            out[cat]["tail_0p5"].append(t5)
            out[cat]["n"].append(counts.sum())
    for c in CAT_NAMES:
        for k in list(out[c].keys()):
            out[c][k] = np.array(out[c][k])
    return out


def make_summary_inclusive(f):
    R_edges = np.linspace(R_MIN, R_MAX, R_NBIN + 1)
    centers = 0.5 * (R_edges[:-1] + R_edges[1:])
    out = {}
    for cat in CAT_NAMES:
        counts, _, _ = read_hist(f, f"R_E_{cat}_inc")
        mu, rms = weighted_moments(counts, centers)
        p16, p50, p84 = weighted_median_percentiles(counts, centers)
        t8, _ = tail_fraction(counts, centers, 0.8)
        t5, _ = tail_fraction(counts, centers, 0.5)
        g = gauss_fit(counts, centers)
        cb = cb_fit(counts, centers)
        out[cat] = {
            "mean": mu, "rms": rms,
            "median": p50, "p16": p16, "p84": p84,
            "tail_0p8": t8, "tail_0p5": t5,
            "gauss_mu": g["mu"] if g else np.nan,
            "gauss_sigma": g["sigma"] if g else np.nan,
            "gauss_mu_err": g["mu_err"] if g else np.nan,
            "gauss_sigma_err": g["sigma_err"] if g else np.nan,
            "cb_mu": cb["mu"] if cb else np.nan,
            "cb_sigma": cb["sigma"] if cb else np.nan,
            "cb_alpha": cb["alpha"] if cb else np.nan,
            "cb_n": cb["n"] if cb else np.nan,
            "cb_mu_err": cb["mu_err"] if cb else np.nan,
            "cb_sigma_err": cb["sigma_err"] if cb else np.nan,
            "n": counts.sum(),
            "counts": counts, "centers": centers,
        }
    return out


def plot_shape_grid(f, summary, name, outpath):
    """Grid of R_E shapes: one panel per pT bin, categories overlaid."""
    n_pt = len(PT_BINS) - 1
    ncol = 4
    nrow = int(np.ceil(n_pt / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.8 * nrow),
                             sharex=True, sharey=False)
    axes = axes.ravel()
    for ip in range(n_pt):
        ax = axes[ip]
        any_plotted = False
        for (cat, label, color, marker, ls) in CATS:
            counts, edges, centers = read_hist(f, f"{name}_{cat}_pt{ip}")
            tot = counts.sum()
            if tot < 200:  # hide low-stat categories in a given bin
                continue
            norm = counts / tot / (edges[1] - edges[0])
            ax.step(centers, norm, where="mid", color=color, linestyle=ls, lw=1.4,
                    label=f"{label} (N={int(tot)})")
            any_plotted = True
        if not any_plotted:
            ax.set_visible(False)
            continue
        ax.set_xlim(0.4, R_MAX)
        ax.set_ylim(1e-3, 30)
        ax.set_yscale("log")
        ax.text(0.04, 0.93,
                f"$p_T^{{truth}} \\in$[{PT_BINS[ip]:.0f}, {PT_BINS[ip+1]:.0f}] GeV",
                transform=ax.transAxes, ha="left", va="top", fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel(r"$R = E^{cluster}/E^{truth}$" if name == "R_E" else r"$R = E_T^{cluster}/p_T^{truth}$")
        if ip % ncol == 0:
            ax.set_ylabel("normalized / bin")
    # Put a single legend on the first populated panel (first pT bin with entries)
    for ax in axes:
        handles, labels_ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels_, fontsize=7, loc="lower center")
            break
    # Hide unused
    for j in range(n_pt, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(f"Response {name} shape, converted vs unconverted photons", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath)
    plt.close(fig)


def plot_mean_vs_pt(summary, outpath):
    pt_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for (cat, label, color, marker, ls) in CATS:
        if cat == "bad":
            continue  # negligible population; suppress from summary plots
        y = summary[cat]["mean"]
        n = summary[cat]["n"]
        err = summary[cat]["rms"] / np.sqrt(np.clip(n, 1, None))
        msk = n > 30
        ax.errorbar(pt_centers[msk], y[msk], yerr=err[msk],
                    fmt=marker, color=color, linestyle=ls, mfc="none" if cat != "unconv" else color,
                    label=label, capsize=2)
    ax.set_xlabel(r"$p_T^{truth}$ [GeV]")
    ax.set_ylabel(r"$\langle R \rangle = \langle E^{cluster}/E^{truth}\rangle$")
    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="lower right")
    ax.set_ylim(0.4, 1.1)
    sphenix_label(ax)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_sigma_vs_pt(summary, outpath):
    """Two panels: Gaussian sigma (core fit) and Crystal Ball sigma (core+tail fit)."""
    pt_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))

    for (cat, label, color, marker, ls) in CATS:
        if cat == "bad":
            continue  # negligible population; CB fit unreliable with <300 events
        n = summary[cat]["n"]
        y_g = summary[cat]["gauss_sigma"]
        y_cb = summary[cat]["cb_sigma"]
        msk_g = (n > 200) & np.isfinite(y_g)
        msk_cb = (n > 500) & np.isfinite(y_cb)
        axes[0].plot(pt_centers[msk_g], y_g[msk_g],
                     marker=marker, color=color, linestyle=ls,
                     mfc="none" if cat != "unconv" else color, label=label)
        axes[1].plot(pt_centers[msk_cb], y_cb[msk_cb],
                     marker=marker, color=color, linestyle=ls,
                     mfc="none" if cat != "unconv" else color, label=label)

    axes[0].set_xlabel(r"$p_T^{truth}$ [GeV]")
    axes[0].set_ylabel(r"Gaussian $\sigma(R)$ (core fit)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=9, loc="lower left")
    axes[0].set_ylim(0.0, 0.12)
    sphenix_label(axes[0], extra=r"Gaussian fit in [peak$-0.10$, peak$+0.10$]")

    axes[1].set_xlabel(r"$p_T^{truth}$ [GeV]")
    axes[1].set_ylabel(r"Crystal Ball $\sigma(R)$")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9, loc="lower left")
    axes[1].set_ylim(0.0, 0.12)
    sphenix_label(axes[1], extra=r"CB fit in $R \in [0.3, 1.2]$")

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_cb_params_vs_pt(summary, outpath):
    """Crystal Ball tail parameters alpha and n vs pT."""
    pt_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))

    for (cat, label, color, marker, ls) in CATS:
        if cat == "bad":
            continue
        n = summary[cat]["n"]
        a = summary[cat]["cb_alpha"]
        nn = summary[cat]["cb_n"]
        msk = (n > 500) & np.isfinite(a) & np.isfinite(nn)
        axes[0].plot(pt_centers[msk], a[msk],
                     marker=marker, color=color, linestyle=ls,
                     mfc="none" if cat != "unconv" else color, label=label)
        axes[1].plot(pt_centers[msk], nn[msk],
                     marker=marker, color=color, linestyle=ls,
                     mfc="none" if cat != "unconv" else color, label=label)

    axes[0].set_xlabel(r"$p_T^{truth}$ [GeV]")
    axes[0].set_ylabel(r"Crystal Ball $\alpha$ (tail transition, $\sigma$ units)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=9, loc="upper right")
    axes[0].set_ylim(0.0, 4.0)
    sphenix_label(axes[0])

    axes[1].set_xlabel(r"$p_T^{truth}$ [GeV]")
    axes[1].set_ylabel(r"Crystal Ball $n$ (tail power-law index)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9, loc="upper right")
    sphenix_label(axes[1])

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_tail_vs_pt(summary, outpath):
    pt_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for (ax, thr_key, title) in [(axes[0], "tail_0p8", "$R<0.8$"),
                                  (axes[1], "tail_0p5", "$R<0.5$")]:
        for (cat, label, color, marker, ls) in CATS:
            if cat == "bad":
                continue
            y = summary[cat][thr_key]
            n = summary[cat]["n"]
            err = np.sqrt(np.clip(y * (1 - y), 0, None) / np.clip(n, 1, None))
            msk = n > 30
            ax.errorbar(pt_centers[msk], y[msk], yerr=err[msk],
                        fmt=marker, color=color, linestyle=ls, mfc="none" if cat != "unconv" else color,
                        label=label, capsize=2)
        ax.set_xlabel(r"$p_T^{truth}$ [GeV]")
        ax.set_ylabel(f"fraction {title}")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc="lower right")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 2.0)
        sphenix_label(ax)
        ax.set_title(f"Tail fraction {title}", fontsize=10)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_angular(f, outpath):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))
    labels = [("dEta", r"$\Delta\eta=\eta^{cluster}-\eta^{truth}$"),
              ("dPhi", r"$\Delta\phi=\phi^{cluster}-\phi^{truth}$"),
              ("dR", r"$\Delta R$")]
    for ax, (name, xlab) in zip(axes, labels):
        for (cat, label, color, marker, ls) in CATS:
            counts, edges, centers = read_hist(f, f"{name}_{cat}")
            tot = counts.sum()
            if tot < 500:  # hide low-stat category (bad-secondary)
                continue
            norm = counts / tot / (edges[1] - edges[0])
            ax.step(centers, norm, where="mid", color=color, linestyle=ls, lw=1.4,
                    label=f"{label} (N={int(tot)})")
        ax.set_xlabel(xlab)
        ax.set_ylabel("normalized / bin")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_yscale("log")
        sphenix_label(ax)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_match_efficiency(f, outpath):
    """Match efficiency per (pT, cat), summed over eta."""
    pt_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    eff_data = {}
    for (cat, label, color, marker, ls) in CATS:
        num, _, _ = read_2d(f, f"num_{cat}")
        den, _, _ = read_2d(f, f"den_{cat}")
        num_pt = num.sum(axis=1)
        den_pt = den.sum(axis=1)
        eff = np.where(den_pt > 0, num_pt / den_pt, 0.0)
        err = np.sqrt(np.clip(eff * (1 - eff), 0, None) / np.clip(den_pt, 1, None))
        msk = den_pt > 30
        ax.errorbar(pt_centers[msk], eff[msk], yerr=err[msk],
                    fmt=marker, color=color, linestyle=ls, mfc="none" if cat != "unconv" else color,
                    label=f"{label} (total={int(den_pt.sum())})", capsize=2)
        eff_data[cat] = {"eff": eff, "num_pt": num_pt, "den_pt": den_pt}
    ax.set_xlabel(r"$p_T^{truth}$ [GeV]")
    ax.set_ylabel("matching efficiency (cluster $E_T>3$ GeV)")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    sphenix_label(ax)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    return eff_data


def plot_eta_dep(f, outpath):
    """Mean R_E and core sigma vs eta per category."""
    R_edges = np.linspace(R_MIN, R_MAX, R_NBIN + 1)
    centers = 0.5 * (R_edges[:-1] + R_edges[1:])
    eta_centers = 0.5 * (ETA_BINS[:-1] + ETA_BINS[1:])
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for (cat, label, color, marker, ls) in CATS:
        if cat == "bad":
            continue
        means = []
        sigmas = []
        ns = []
        for ie in range(len(ETA_BINS) - 1):
            counts, _, _ = read_hist(f, f"R_E_{cat}_eta{ie}")
            mu, _ = weighted_moments(counts, centers)
            fit = gauss_fit(counts, centers)
            means.append(mu)
            sigmas.append(fit["sigma"] if fit else np.nan)
            ns.append(counts.sum())
        means = np.array(means)
        sigmas = np.array(sigmas)
        ns = np.array(ns)
        msk = ns > 100
        axes[0].plot(eta_centers[msk], means[msk], marker=marker, color=color,
                     linestyle=ls, mfc="none" if cat != "unconv" else color, label=label)
        axes[1].plot(eta_centers[msk], sigmas[msk], marker=marker, color=color,
                     linestyle=ls, mfc="none" if cat != "unconv" else color, label=label)
    axes[0].set_xlabel(r"$\eta^{truth}$")
    axes[0].set_ylabel(r"$\langle R \rangle$")
    axes[0].axhline(1.0, color="gray", linestyle=":", lw=0.8)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=9)
    axes[0].set_ylim(0.5, 1.05)
    sphenix_label(axes[0])
    axes[1].set_xlabel(r"$\eta^{truth}$")
    axes[1].set_ylabel(r"Gaussian $\sigma(R)$ (core fit)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9, loc="lower right")
    axes[1].set_ylim(0.0, 0.12)
    sphenix_label(axes[1])
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def main():
    f = uproot.open(ROOT_PATH)
    summary = make_summary(f)
    inclusive = make_summary_inclusive(f)
    plot_shape_grid(f, summary, "R_E", FIG / "response_shape_grid.pdf")
    plot_shape_grid(f, summary, "R_ET", FIG / "response_shape_ET_grid.pdf")
    plot_mean_vs_pt(summary, FIG / "response_mean_vs_pt.pdf")
    plot_sigma_vs_pt(summary, FIG / "response_sigma_vs_pt.pdf")
    plot_cb_params_vs_pt(summary, FIG / "response_cb_params_vs_pt.pdf")
    plot_tail_vs_pt(summary, FIG / "response_tailfrac_vs_pt.pdf")
    plot_angular(f, FIG / "response_angular_deltas.pdf")
    eff_data = plot_match_efficiency(f, FIG / "response_match_eff_vs_pt.pdf")
    plot_eta_dep(f, FIG / "response_eta_dep.pdf")

    # Save numeric summary
    out = {
        "pt_edges": PT_BINS,
        "eta_edges": ETA_BINS,
        "summary_per_pt": {c: {k: v.tolist() for k, v in summary[c].items()} for c in CAT_NAMES},
        "summary_inclusive": {c: {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                                   for k, v in inclusive[c].items() if k not in ("counts", "centers")}
                              for c in CAT_NAMES},
        "match_eff": {c: {"eff_vs_pt": eff_data[c]["eff"].tolist(),
                          "n_truth_vs_pt": eff_data[c]["den_pt"].tolist(),
                          "n_matched_vs_pt": eff_data[c]["num_pt"].tolist()}
                      for c in CAT_NAMES},
    }
    with open(HERE / "rootFiles" / "summary.pkl", "wb") as sf:
        pickle.dump(out, sf)
    # Also dump a text summary
    with open(HERE / "rootFiles" / "summary.txt", "w") as tf:
        tf.write("Inclusive (pT 10-34, |eta|<0.7) response summary\n")
        tf.write("=" * 60 + "\n")
        for c in CAT_NAMES:
            s = inclusive[c]
            tf.write(f"\n{c}:\n")
            tf.write(f"  N matched = {int(s['n'])}\n")
            tf.write(f"  <R_E>  = {s['mean']:.4f}\n")
            tf.write(f"  RMS    = {s['rms']:.4f}\n")
            tf.write(f"  median = {s['median']:.4f}  [p16={s['p16']:.4f}, p84={s['p84']:.4f}]\n")
            tf.write(f"  Gauss (+-0.10 core) mu = {s['gauss_mu']:.4f}, sigma = {s['gauss_sigma']:.4f}\n")
            tf.write(f"  Crystal Ball        mu = {s['cb_mu']:.4f}, sigma = {s['cb_sigma']:.4f}, "
                     f"alpha = {s['cb_alpha']:.3f}, n = {s['cb_n']:.3f}\n")
            tf.write(f"  frac R<0.8 = {s['tail_0p8']:.4f}\n")
            tf.write(f"  frac R<0.5 = {s['tail_0p5']:.4f}\n")
        tf.write("\nMatch efficiency (integrated 10-34 GeV):\n")
        for c in CAT_NAMES:
            tot_num = np.array(out["match_eff"][c]["n_matched_vs_pt"]).sum()
            tot_den = np.array(out["match_eff"][c]["n_truth_vs_pt"]).sum()
            eff = tot_num / tot_den if tot_den > 0 else np.nan
            tf.write(f"  {c}: {tot_num}/{tot_den} = {eff:.4f}\n")
    print("[make_plots] wrote figures to", FIG)
    print("[make_plots] summary:")
    for c in CAT_NAMES:
        s = inclusive[c]
        print(f"  {c}: <R>={s['mean']:.4f}, gauss mu={s['gauss_mu']:.4f}, "
              f"gauss sigma={s['gauss_sigma']:.4f}, "
              f"CB mu={s['cb_mu']:.4f}, CB sigma={s['cb_sigma']:.4f}, "
              f"CB alpha={s['cb_alpha']:.3f}, CB n={s['cb_n']:.3f}, "
              f"tail(R<0.8)={s['tail_0p8']:.4f}, N={int(s['n'])}")


if __name__ == "__main__":
    main()
