"""Alt check 1: regrid tower maps into 2x2 / 4x4 super-bins, redo log(R) fit
and mask-count analysis. Writes comparison table + 3-panel plot (1x1 / 2x2 /
4x4) of log(R) distributions per level.
"""
import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.stats import norm

RES = '/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results'
FIG = '/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures'

LEVELS = ['preselect', 'common', 'tight', 'tight_iso']
GROUPS = [1, 2, 4]

def fit_gauss_3sigma(x):
    m, s = x.mean(), x.std()
    for _ in range(3):
        k = np.abs(x - m) < 3 * s
        if k.sum() == 0: break
        m, s = x[k].mean(), x[k].std()
    return m, s

def rebin2d(a, n):
    nx, ny = a.shape
    assert nx % n == 0 and ny % n == 0
    return a.reshape(nx // n, n, ny // n, n).sum(axis=(1, 3))

f = uproot.open(f'{RES}/Photon_final_bdt_nom.root')
print(f"{'level':<12s} {'grp':>4s} {'nbins':>6s} {'mu_logR':>8s} {'sig_logR':>10s} {'hd':>4s} {'zlt2':>6s} {'zlt5':>6s} {'def_zlt2_%':>10s}")
for lvl in LEVELS:
    hmc_full = f[f'h_etaphi_tower_{lvl}_mc_inclusive'].values()
    hda_full = f[f'h_etaphi_tower_{lvl}_data'].values()
    for n in GROUPS:
        hmc = rebin2d(hmc_full, n)
        hda = rebin2d(hda_full, n)
        Nm, Nd = hmc.sum(), hda.sum()
        # mask of bins with both nonzero
        both = (hmc > 0) & (hda > 0)
        R = np.where(both, (hda / Nd) / (hmc / Nm), 0.0)
        logR_valid = np.log(R[both])
        mu, sig = fit_gauss_3sigma(logR_valid)
        # categories
        z = np.full_like(R, np.nan)
        z[both] = (np.log(R[both]) - mu) / sig
        n_hd   = int((hmc > 0).sum() - both.sum())  # MC>0, data=0
        n_zlt5 = int((z < -5).sum())
        n_zlt2 = int(((z < -2) & (z >= -5)).sum())
        # deficit z<-2
        mu_exp = hmc * (Nd / Nm)
        deficit = np.where((z < -2), np.maximum(0, mu_exp - hda), 0).sum()
        deficit_pct = 100 * deficit / Nd
        print(f"{lvl:<12s} {n:>3d}x{n} {R.size:>6d} {mu:>8.3f} {sig:>10.3f} {n_hd:>4d} {n_zlt2:>6d} {n_zlt5:>6d} {deficit_pct:>10.2f}")

# -- plot 3-panel per level ------------------------------------------------
fig, axes = plt.subplots(4, 3, figsize=(13, 14), sharex=False, sharey=False)
for il, lvl in enumerate(LEVELS):
    hmc_full = f[f'h_etaphi_tower_{lvl}_mc_inclusive'].values()
    hda_full = f[f'h_etaphi_tower_{lvl}_data'].values()
    for ig, n in enumerate(GROUPS):
        ax = axes[il, ig]
        hmc = rebin2d(hmc_full, n); hda = rebin2d(hda_full, n)
        Nm, Nd = hmc.sum(), hda.sum()
        both = (hmc > 0) & (hda > 0)
        R = (hda[both] / Nd) / (hmc[both] / Nm)
        lR = np.log(R[R > 0])
        mu, sig = fit_gauss_3sigma(lR)
        bins = np.linspace(-4, 4, 161)
        ax.hist(lR, bins=bins, histtype='step', color='k', lw=1.5)
        # gauss overlay
        x = np.linspace(-4, 4, 400)
        bw = bins[1] - bins[0]
        gv = len(lR) * bw * norm.pdf(x, mu, sig)
        ax.plot(x, gv, 'r--', lw=1.5, label=f'Gauss fit\n$\\mu$={mu:.2f}, $\\sigma$={sig:.2f}')
        ax.axvline(mu, color='r', ls=':', lw=1)
        ax.set_yscale('log')
        ax.set_ylim(0.8, None)
        ax.set_xlim(-4, 4)
        ax.set_xlabel('log(R)' if il == 3 else '')
        ax.set_ylabel(f'{lvl}\nbins' if ig == 0 else '')
        ax.set_title(f'{n}x{n} tower grouping' if il == 0 else '')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(alpha=0.3)
fig.suptitle('log(R) distributions at 1x1 / 2x2 / 4x4 tower grouping', y=0.995, fontsize=12)
fig.tight_layout()
fig.savefig(f'{FIG}/alt_check_rebin.pdf', bbox_inches='tight')
fig.savefig(f'{FIG}/alt_check_rebin.png', dpi=150, bbox_inches='tight')
print(f'\nwrote {FIG}/alt_check_rebin.pdf')
