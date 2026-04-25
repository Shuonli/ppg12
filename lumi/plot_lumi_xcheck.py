#!/usr/bin/env python3
"""
plot_lumi_xcheck.py

Read lumi/lumi_from_slimtree_allz.list + lumi/trig29_30_31_allz.list and
produce:
  (a) per-run scatter L_UC_slim vs L_UC_joey
  (b) per-run scatter L_Corr_slim vs L_Corr_joey  (F_corr taken from Joey)
  (c) per-run ratio slim/joey vs run number (UC and Corr panels)
  (d) summary text: aggregate totals, crossing-angle split, outlier runs.

Output figure: plotting/figures/lumi_slimtree_xcheck.pdf
Output summary: reports/lumi_slimtree_xcheck_summary.txt
"""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


RUN_XING_SPLIT = 51274  # 0 mrad below, 1.5 mrad at-or-above (nominal)


def load_slim(path):
    """Return dict rn -> (L_UC_slim, L_Corr_slim, L_Corr_joey, L_UC_joey,
    PS_phot, max_live_mbd, min_live_mbd, nevt)."""
    out = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            try:
                rn = int(parts[0])
            except ValueError:
                continue
            L_uc_slim = float(parts[1])
            L_corr_slim = float(parts[2])
            L_corr_joey = float(parts[3])
            L_uc_joey = float(parts[4])
            ps_phot = float(parts[5])
            max_live_mbd = int(parts[7])
            min_live_mbd = int(parts[8])
            nevt = int(parts[-1])
            out[rn] = dict(L_uc_slim=L_uc_slim, L_corr_slim=L_corr_slim,
                            L_corr_joey=L_corr_joey, L_uc_joey=L_uc_joey,
                            ps_phot=ps_phot, max_live_mbd=max_live_mbd,
                            min_live_mbd=min_live_mbd, nevt=nevt)
    return out


def summarize(slim, fout):
    rns = sorted(slim)
    def tot(key, subset=None):
        s = 0.0
        for rn in (subset or rns):
            s += slim[rn][key]
        return s

    rns_0mrad = [rn for rn in rns if rn < RUN_XING_SPLIT]
    rns_15mrad = [rn for rn in rns if rn >= RUN_XING_SPLIT]

    lines = []
    lines.append('slimtree-based all-z luminosity cross-check (Photon_4_GeV = bit 30)')
    lines.append('=' * 72)
    lines.append(f'Runs in slimtree output: {len(rns)}')
    lines.append(f'  0 mrad   (RN <  {RUN_XING_SPLIT}): {len(rns_0mrad)} runs')
    lines.append(f'  1.5 mrad (RN >= {RUN_XING_SPLIT}): {len(rns_15mrad)} runs')
    lines.append('')
    lines.append('Aggregate totals (pb^-1):')
    for label, subset in [('All', rns), ('0 mrad', rns_0mrad), ('1.5 mrad', rns_15mrad)]:
        Luc_s = tot('L_uc_slim', subset)
        Lcorr_s = tot('L_corr_slim', subset)
        Luc_j = tot('L_uc_joey', subset)
        Lcorr_j = tot('L_corr_joey', subset)
        lines.append(f'  {label:9s}: UC   slim={Luc_s:9.4f}  joey={Luc_j:9.4f}  ratio={Luc_s/Luc_j if Luc_j else float("nan"):.4f}')
        lines.append(f'  {label:9s}: Corr slim={Lcorr_s:9.4f}  joey={Lcorr_j:9.4f}  ratio={Lcorr_s/Lcorr_j if Lcorr_j else float("nan"):.4f}')
    lines.append('')

    uc_ratios = []
    corr_ratios = []
    for rn in rns:
        if slim[rn]['L_uc_joey'] > 0:
            uc_ratios.append(slim[rn]['L_uc_slim'] / slim[rn]['L_uc_joey'])
        if slim[rn]['L_corr_joey'] > 0:
            corr_ratios.append(slim[rn]['L_corr_slim'] / slim[rn]['L_corr_joey'])
    uc_ratios = np.array(uc_ratios)
    corr_ratios = np.array(corr_ratios)
    lines.append(f'Per-run ratio slim/joey (UC):   median={np.median(uc_ratios):.4f}  mean={np.mean(uc_ratios):.4f}  std={np.std(uc_ratios):.4f}')
    lines.append(f'Per-run ratio slim/joey (Corr): median={np.median(corr_ratios):.4f}  mean={np.mean(corr_ratios):.4f}  std={np.std(corr_ratios):.4f}')
    lines.append(f'Runs with UC ratio in [0.95, 1.05]: {int(((uc_ratios >= 0.95) & (uc_ratios <= 1.05)).sum())} / {len(uc_ratios)}')
    lines.append(f'Runs with UC ratio > 1.10:          {int((uc_ratios > 1.10).sum())} / {len(uc_ratios)}')
    lines.append(f'Runs with UC ratio < 0.90:          {int((uc_ratios < 0.90).sum())} / {len(uc_ratios)}')
    lines.append('')

    # Worst outliers
    rn_arr = np.array([rn for rn in rns if slim[rn]['L_uc_joey'] > 0])
    ratio_arr = np.array([slim[rn]['L_uc_slim']/slim[rn]['L_uc_joey'] for rn in rn_arr])
    lines.append('Top 10 UC outliers (|ratio - 1| largest):')
    lines.append('  RN      L_UC_slim    L_UC_joey    ratio    PS_phot  nevt_slim')
    order = np.argsort(-np.abs(ratio_arr - 1.0))
    for idx in order[:10]:
        rn = int(rn_arr[idx])
        s = slim[rn]
        lines.append(f'  {rn}  {s["L_uc_slim"]:10.6f}  {s["L_uc_joey"]:10.6f}  {ratio_arr[idx]:7.3f}  {s["ps_phot"]:7.2f}  {s["nevt"]}')

    text = '\n'.join(lines) + '\n'
    if fout:
        with open(fout, 'w') as f:
            f.write(text)
    print(text)
    return uc_ratios, corr_ratios


def make_plots(slim, pdf_path):
    os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
    rns = sorted(slim)
    L_uc_slim = np.array([slim[rn]['L_uc_slim'] for rn in rns])
    L_uc_joey = np.array([slim[rn]['L_uc_joey'] for rn in rns])
    L_corr_slim = np.array([slim[rn]['L_corr_slim'] for rn in rns])
    L_corr_joey = np.array([slim[rn]['L_corr_joey'] for rn in rns])
    rn_arr = np.array(rns)
    color_xing = np.where(rn_arr < RUN_XING_SPLIT, 'C0', 'C3')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) UC scatter
    ax = axes[0, 0]
    ax.scatter(L_uc_joey, L_uc_slim, c=color_xing, s=15, alpha=0.7)
    mx = max(L_uc_joey.max(), L_uc_slim.max()) * 1.05
    ax.plot([0, mx], [0, mx], 'k--', lw=1, label='y=x')
    ax.set_xlabel(r'$L_{UC}$ from Joey (Bit30UC)  [pb$^{-1}$]')
    ax.set_ylabel(r'$L_{UC}$ from slimtree  [pb$^{-1}$]')
    ax.set_title(f'all-z UC lumi per run (sum slim={L_uc_slim.sum():.3f}, joey={L_uc_joey.sum():.3f} pb$^{{-1}}$)')
    ax.legend(loc='upper left')
    ax.grid(alpha=0.3)

    # (b) Corr scatter
    ax = axes[0, 1]
    ax.scatter(L_corr_joey, L_corr_slim, c=color_xing, s=15, alpha=0.7)
    mx = max(L_corr_joey.max(), L_corr_slim.max()) * 1.05
    ax.plot([0, mx], [0, mx], 'k--', lw=1, label='y=x')
    ax.set_xlabel(r'$L_{Corr}$ from Joey (Bit30Corr)  [pb$^{-1}$]')
    ax.set_ylabel(r'$L_{Corr}$ from slimtree  [pb$^{-1}$]  (= $L_{UC}^{slim} \cdot F_{corr}^{joey}$)')
    ax.set_title(f'all-z Corr lumi per run (sum slim={L_corr_slim.sum():.3f}, joey={L_corr_joey.sum():.3f} pb$^{{-1}}$)')
    ax.legend(loc='upper left')
    ax.grid(alpha=0.3)

    # (c) ratio vs RN, UC
    ax = axes[1, 0]
    mask = L_uc_joey > 0
    ratio_uc = np.full_like(L_uc_slim, np.nan)
    ratio_uc[mask] = L_uc_slim[mask] / L_uc_joey[mask]
    ax.scatter(rn_arr, ratio_uc, c=color_xing, s=15, alpha=0.7)
    ax.axhline(1.0, color='k', ls='--', lw=1)
    ax.axvline(RUN_XING_SPLIT, color='gray', ls=':', lw=1, label=f'xing split RN={RUN_XING_SPLIT}')
    ax.set_xlabel('Run number')
    ax.set_ylabel(r'$L_{UC}^{slim} / L_{UC}^{joey}$')
    ax.set_title('Per-run UC ratio')
    ax.set_ylim(0, 3)
    ax.legend(loc='upper left')
    ax.grid(alpha=0.3)

    # (d) ratio vs RN, Corr
    ax = axes[1, 1]
    mask = L_corr_joey > 0
    ratio_corr = np.full_like(L_corr_slim, np.nan)
    ratio_corr[mask] = L_corr_slim[mask] / L_corr_joey[mask]
    ax.scatter(rn_arr, ratio_corr, c=color_xing, s=15, alpha=0.7)
    ax.axhline(1.0, color='k', ls='--', lw=1)
    ax.axvline(RUN_XING_SPLIT, color='gray', ls=':', lw=1, label=f'xing split RN={RUN_XING_SPLIT}')
    ax.set_xlabel('Run number')
    ax.set_ylabel(r'$L_{Corr}^{slim} / L_{Corr}^{joey}$')
    ax.set_title('Per-run Corr ratio')
    ax.set_ylim(0, 3)
    ax.legend(loc='upper left')
    ax.grid(alpha=0.3)

    fig.suptitle('slimtree-based all-z luminosity cross-check (Photon_4_GeV, $\sigma_{MBD}=25.2$ mb)',
                 fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(pdf_path, bbox_inches='tight')
    print(f'Wrote {pdf_path}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slim-list',
                    default='/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/lumi/lumi_from_slimtree_allz.list')
    ap.add_argument('--pdf',
                    default='/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/lumi_slimtree_xcheck.pdf')
    ap.add_argument('--summary',
                    default='/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/lumi_slimtree_xcheck_summary.txt')
    args = ap.parse_args()

    slim = load_slim(args.slim_list)
    if not slim:
        raise SystemExit(f'No rows in {args.slim_list}')
    os.makedirs(os.path.dirname(args.summary), exist_ok=True)
    summarize(slim, args.summary)
    make_plots(slim, args.pdf)


if __name__ == '__main__':
    main()
