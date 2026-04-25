"""Dump detailed per-level statistics from the official tower-map pipeline.

Input : efficiencytool/results/Photon_final_bdt_nom.root
Output: stdout tables
"""
import numpy as np
import uproot

F = '/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root'
LEVELS = ['preselect', 'common', 'tight', 'tight_iso']

# fiducial in tower index: |eta|<0.7 <=> ieta in [17, 78]
IETA_FID_LO = 17
IETA_FID_HI = 78

f = uproot.open(F)

def read(name):
    h = f[name]
    counts = h.values()
    try:
        variances = h.variances()
    except Exception:
        variances = counts
    return counts, variances

def fiducial_mask(nx=96, ny=256):
    m = np.zeros((nx, ny), dtype=bool)
    m[IETA_FID_LO:IETA_FID_HI+1, :] = True
    return m

def analyse(lvl, restrict_fiducial=False):
    n_mc, v_mc = read(f'h_etaphi_tower_{lvl}_mc_inclusive')
    n_da, v_da = read(f'h_etaphi_tower_{lvl}_data')
    mask = fiducial_mask() if restrict_fiducial else np.ones_like(n_mc, dtype=bool)
    n_mc = n_mc * mask; n_da = n_da * mask
    v_mc = v_mc * mask

    N_mc = n_mc.sum(); N_da = n_da.sum()
    # effective N per bin for MC (weight^2 based)
    with np.errstate(divide='ignore', invalid='ignore'):
        neff_mc = np.where(v_mc > 0, n_mc**2 / v_mc, 0.0)

    tag = 'FIDUCIAL (|eta|<0.7)' if restrict_fiducial else 'FULL RANGE'
    print(f'\n=== {lvl.upper()} ({tag}) ===')
    print(f'  N_mc_inclusive (weighted) = {N_mc:.4e}')
    print(f'  N_data                    = {N_da:.4e}')
    # bin counts
    nonzero_mc = int((n_mc > 0).sum())
    nonzero_da = int((n_da > 0).sum())
    zero_mc_zero_da = int(((n_mc <= 0) & (n_da <= 0) & mask).sum())
    nonzero_mc_zero_da = int(((n_mc > 0) & (n_da <= 0)).sum())
    zero_mc_nonzero_da = int(((n_mc <= 0) & (n_da > 0)).sum())
    fiducial_bins = int(mask.sum())
    print(f'  fiducial bins      = {fiducial_bins}')
    print(f'  MC nonzero bins    = {nonzero_mc}  ({100*nonzero_mc/fiducial_bins:.2f}%)')
    print(f'  data nonzero bins  = {nonzero_da}  ({100*nonzero_da/fiducial_bins:.2f}%)')
    print(f'  both empty         = {zero_mc_zero_da}')
    print(f'  MC live, data=0    = {nonzero_mc_zero_da}  <-- HARD DEAD IN DATA')
    print(f'  MC=0, data live    = {zero_mc_nonzero_da}')

    # z-score only where nm>0
    p_mc = np.where(N_mc > 0, n_mc / max(N_mc, 1e-30), 0.0)
    p_da = np.where(N_da > 0, n_da / max(N_da, 1e-30), 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        R = np.where(p_mc > 0, p_da / p_mc, np.nan)
        relerr = np.sqrt(1.0 / np.maximum(n_da, 1.0) + 1.0 / np.maximum(neff_mc, 1.0))
        z = np.where((p_mc > 0) & np.isfinite(R), (R - 1) / (R * relerr), np.nan)

    # MC effective N statistics (for bins with MC content)
    neff_arr = neff_mc[n_mc > 0]
    print(f'  MC Neff per bin    : median={np.median(neff_arr):.1f}  p10={np.percentile(neff_arr,10):.1f}  p90={np.percentile(neff_arr,90):.1f}')

    # z-score tail distribution
    zvalid = z[np.isfinite(z)]
    n_z = len(zvalid)
    for thr in [-2, -3, -5, -10, 2, 3, 5]:
        if thr < 0:
            count = int((zvalid < thr).sum())
            side = f'z<{thr}'
        else:
            count = int((zvalid > thr).sum())
            side = f'z>{thr}'
        print(f'  {side:7s}: {count:5d}  ({100*count/n_z:.3f} % of valid bins)')

    # Deficit accounting (MC-normalized-to-data scale)
    mu = n_mc * (N_da / max(N_mc, 1e-30))
    # z<-2 bins
    mask_z2 = (z < -2) & np.isfinite(z)
    mask_z5 = (z < -5) & np.isfinite(z)
    def_z2 = np.sum(np.maximum(0.0, mu - n_da) * mask_z2)
    def_z5 = np.sum(np.maximum(0.0, mu - n_da) * mask_z5)
    excess_zp2 = np.sum(np.maximum(0.0, n_da - mu) * ((z > 2) & np.isfinite(z)))
    print(f'  data deficit z<-2 = {def_z2:.1f} clusters  ({100*def_z2/N_da:.2f} % of data)')
    print(f'  data deficit z<-5 = {def_z5:.1f} clusters  ({100*def_z5/N_da:.2f} % of data)')
    print(f'  data excess  z>+2 = {excess_zp2:.1f} clusters  ({100*excess_zp2/N_da:.2f} % of data)')
    print(f'  net: (def_z2 - excess_zp2)/N_data = {100*(def_z2 - excess_zp2)/N_da:+.2f} %')

    # Top-5 most deficient bins (by absolute mu - nd) — where are the holes?
    deficit_arr = np.where(mask_z2, mu - n_da, 0.0)
    flat_idx = np.argsort(deficit_arr.ravel())[::-1][:5]
    print('  Top-5 dead-in-data towers (z<-2):')
    print(f'  {"ieta":>6} {"iphi":>6} {"n_data":>10} {"n_MC":>12} {"mu":>10} {"z":>8}')
    for idx in flat_idx:
        ix, iy = np.unravel_index(idx, n_mc.shape)
        if deficit_arr[ix, iy] <= 0:
            break
        print(f'  {ix:>6} {iy:>6} {n_da[ix,iy]:>10.1f} {n_mc[ix,iy]:>12.3e} {mu[ix,iy]:>10.1f} {z[ix,iy]:>8.1f}')

    # Top-5 z>+2 bins (data excess)
    excess_arr = np.where((z > 2) & np.isfinite(z), n_da - mu, 0.0)
    flat_idx = np.argsort(excess_arr.ravel())[::-1][:5]
    print('  Top-5 data-excess towers (z>+2):')
    print(f'  {"ieta":>6} {"iphi":>6} {"n_data":>10} {"n_MC":>12} {"mu":>10} {"z":>8} {"Neff":>8}')
    for idx in flat_idx:
        ix, iy = np.unravel_index(idx, n_mc.shape)
        if excess_arr[ix, iy] <= 0:
            break
        print(f'  {ix:>6} {iy:>6} {n_da[ix,iy]:>10.1f} {n_mc[ix,iy]:>12.3e} {mu[ix,iy]:>10.1f} {z[ix,iy]:>8.1f} {neff_mc[ix,iy]:>8.1f}')


for lvl in LEVELS:
    for restrict in (False, True):
        analyse(lvl, restrict_fiducial=restrict)
