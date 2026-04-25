"""Alt check 2: data-only phi-symmetry. For each eta row in the tower map,
assume data is flat in phi (physics is phi-symmetric at fixed eta). Fit a
Gaussian to the distribution of counts across phi; flag phi bins with
z<-2 (count far below mean) as dead candidates. Compare to MC-based masks.
"""
import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import norm

RES = '/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results'
FIG = '/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures'

LEVELS = ['preselect', 'common', 'tight', 'tight_iso']

def fit_gauss_3sigma(x):
    m, s = x.mean(), x.std()
    for _ in range(3):
        k = np.abs(x - m) < 3 * s
        if k.sum() == 0: break
        m, s = x[k].mean(), x[k].std()
    return m, s

# Per-eta band grouping: 2 eta rows combined -> 48 bands of 2 eta x 256 phi
# (boosts per-band stats vs row-by-row, which is too sparse at tight)
ETA_GROUP = 2  # group 2 eta rows
# Fiducial ieta cut: bands fully inside [17, 78] require band >= 9 & <= 38
# (since band index b spans eta rows [2b, 2b+1])
IETA_FID_LO, IETA_FID_HI = 17, 78

f = uproot.open(f'{RES}/Photon_final_bdt_nom.root')

print(f"{'level':<12s} {'n_dead_z<-2':>14s} {'n_dead_z<-5':>14s} "
      f"{'hd_for_ref':>12s} {'deficit_%':>10s}")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for il, lvl in enumerate(LEVELS):
    hda_full = f[f'h_etaphi_tower_{lvl}_data'].values()
    # For reference: hard-dead count (data=0 where MC>0)
    hmc_full = f[f'h_etaphi_tower_{lvl}_mc_inclusive'].values()

    neta, nphi = hda_full.shape  # 96 x 256
    # Group eta rows by ETA_GROUP (2)
    n_bands = neta // ETA_GROUP
    banded = hda_full.reshape(n_bands, ETA_GROUP, nphi).sum(axis=1)

    # For each band, compute per-phi counts; fit Gaussian to those counts;
    # flag phi bins with z<-2 (below mean by 2 sigma) within that band.
    dead_mask = np.zeros((n_bands, nphi), dtype=bool)
    deadv_mask = np.zeros((n_bands, nphi), dtype=bool)
    total_data = 0
    total_deficit_zlt2 = 0
    total_deficit_zlt5 = 0
    for b in range(n_bands):
        row = banded[b].astype(float)
        if row.sum() < 10:
            continue  # too few events; skip
        mu, sig = fit_gauss_3sigma(row)
        if sig <= 0:
            continue
        z = (row - mu) / sig
        dead_mask[b, z < -2] = True
        deadv_mask[b, z < -5] = True
        # deficit = (mu - row) for z<-2 bins (positive)
        deficit = np.where(z < -2, np.maximum(0, mu - row), 0).sum()
        deficit_v = np.where(z < -5, np.maximum(0, mu - row), 0).sum()
        total_data += row.sum()
        total_deficit_zlt2 += deficit
        total_deficit_zlt5 += deficit_v

    # Expand back to 96x256 to report per-tower counts for fiducial
    dead_tower = np.repeat(dead_mask, ETA_GROUP, axis=0)
    deadv_tower = np.repeat(deadv_mask, ETA_GROUP, axis=0)
    # Fiducial
    fid = np.zeros((neta, nphi), dtype=bool)
    fid[IETA_FID_LO:IETA_FID_HI + 1, :] = True
    n_dead_fid = int((dead_tower & fid).sum())
    n_deadv_fid = int((deadv_tower & fid).sum())
    # Hard-dead reference
    hd = ((hmc_full > 0) & (hda_full == 0))
    n_hd = int((hd & fid).sum())

    deficit_pct = 100 * total_deficit_zlt2 / total_data if total_data > 0 else 0
    print(f"{lvl:<12s} {n_dead_fid:>14d} {n_deadv_fid:>14d} {n_hd:>12d} {deficit_pct:>10.2f}")

    # Plot: show data count map + overlay dead_tower mask as red 'X' markers
    ax = axes[il]
    im = ax.imshow(hda_full.T, origin='lower', aspect='auto',
                   extent=(0, 96, 0, 256), cmap='viridis',
                   norm=plt.matplotlib.colors.LogNorm(vmin=max(1, hda_full.min()),
                                                     vmax=hda_full.max()))
    # Mark dead phi bins (from 2x1 bands -> repeat into towers)
    where_dead = np.where(dead_tower)
    ax.scatter(where_dead[0] + 0.5, where_dead[1] + 0.5,
               s=1.5, c='red', marker='s', alpha=0.6, label='z<-2 (phi-symm)')
    # Fiducial rectangle
    ax.add_patch(Rectangle((IETA_FID_LO, 0), IETA_FID_HI - IETA_FID_LO + 1, nphi,
                            fill=False, edgecolor='cyan', lw=2, ls='--'))
    ax.set_xlabel('cluster i$\\eta$')
    ax.set_ylabel('cluster i$\\phi$')
    ax.set_title(f'{lvl}: z<-2 (phi-symm, data-only) = {n_dead_fid} tower-cells\n'
                 f'(hard-dead ref = {n_hd}, deficit = {deficit_pct:.2f}% of data)')
    ax.legend(loc='upper right', fontsize=8)

fig.suptitle('Alt check 2: data-only phi-symmetry dead-tower finder\n'
             '(group 2 eta rows; fit Gaussian to per-phi counts; flag z<-2)',
             fontsize=12)
fig.tight_layout()
fig.savefig(f'{FIG}/alt_check_phisymm.pdf', bbox_inches='tight')
fig.savefig(f'{FIG}/alt_check_phisymm.png', dpi=150, bbox_inches='tight')
print(f'\nwrote {FIG}/alt_check_phisymm.pdf')
