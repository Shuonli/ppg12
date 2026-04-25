#!/usr/bin/env python3
"""
lumi_from_slimtree.py

Compute per-run all-z integrated luminosity from the analysis slimtree's
GL1 scaler branches, with and without pileup correction.

UC (uncorrected):
    L_UC(run) = max_i{currentscaler_live[MBD_BIT]}  /  PS_phot[run]  /  (sigma_MBD * 1e9)
where MBD_BIT = 10 (MBD_NandS_geq_1), PS_phot = trigger_prescale[PHOT_BIT],
sigma_MBD = 25.2 mb, and the 1e9 factor converts mb -> pb.

Corr (pileup-corrected):
    L_Corr(run) = L_UC(run) * F_corr(run)
with F_corr pulled per run from Joey's list as F_corr = Bit30Corr / Bit30UC.

Only events written to the slimtree contribute. Because CaloAna24 only
fills slimtree entries for events with >=1 cluster (data additionally has
|vertexz|<60 cm), max(currentscaler_live[MBD_BIT]) is the MBD scaler
value at the last selected event, which is typically very close to
end-of-run. Min is also tracked for QA and optional "max-min" cross
checks.
"""

import argparse
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import uproot


SIGMA_MBD_MB = 25.2
PB_PER_MB = 1e9

MBD_BIT = 10
PHOT_BIT = 30

INIT_STATS = dict(max_live_mbd=-1, min_live_mbd=-1,
                  max_raw_mbd=-1, max_scaled_mbd=-1,
                  max_live_phot=-1, ps_phot=-1.0, ps_mbd=-1.0,
                  nevt=0)


def blank_stats():
    return dict(INIT_STATS)


def update_stats(s, field_max, field_min, val):
    if field_max is not None:
        if s[field_max] < 0 or val > s[field_max]:
            s[field_max] = int(val)
    if field_min is not None:
        if s[field_min] < 0 or val < s[field_min]:
            s[field_min] = int(val)


def per_file_scan(fpath):
    stats = {}
    try:
        with uproot.open(fpath) as f:
            t = f['slimtree']
            for arrays in t.iterate(
                ['runnumber', 'currentscaler_live', 'currentscaler_raw',
                 'currentscaler_scaled', 'trigger_prescale'],
                library='np', step_size=200_000,
            ):
                rn = arrays['runnumber']
                live = arrays['currentscaler_live']
                raw = arrays['currentscaler_raw']
                scaled = arrays['currentscaler_scaled']
                ps = arrays['trigger_prescale']
                uniq = np.unique(rn)
                for rn_val in uniq:
                    rn_int = int(rn_val)
                    mask = (rn == rn_val)
                    s = stats.setdefault(rn_int, blank_stats())
                    lmbd = live[mask, MBD_BIT]
                    rmbd = raw[mask, MBD_BIT]
                    smbd = scaled[mask, MBD_BIT]
                    lphot = live[mask, PHOT_BIT]
                    update_stats(s, 'max_live_mbd', 'min_live_mbd', lmbd.max())
                    if lmbd.size:
                        update_stats(s, None, 'min_live_mbd', lmbd.min())
                    update_stats(s, 'max_raw_mbd', None, rmbd.max())
                    update_stats(s, 'max_scaled_mbd', None, smbd.max())
                    update_stats(s, 'max_live_phot', None, lphot.max())
                    s['nevt'] += int(mask.sum())
                    if s['ps_phot'] <= 0:
                        v = ps[mask, PHOT_BIT]
                        v = v[np.isfinite(v) & (v > 0)]
                        if v.size:
                            s['ps_phot'] = float(v[0])
                    if s['ps_mbd'] <= 0:
                        v = ps[mask, MBD_BIT]
                        v = v[np.isfinite(v) & (v > 0)]
                        if v.size:
                            s['ps_mbd'] = float(v[0])
    except Exception as e:
        print(f'ERROR reading {fpath}: {e}', file=sys.stderr)
    return stats


def merge_stats(all_stats):
    merged = {}
    for stats in all_stats:
        for rn, s in stats.items():
            m = merged.setdefault(rn, blank_stats())
            for k in ('max_live_mbd', 'max_raw_mbd', 'max_scaled_mbd', 'max_live_phot'):
                if s[k] > m[k]:
                    m[k] = s[k]
            if s['min_live_mbd'] >= 0:
                if m['min_live_mbd'] < 0 or s['min_live_mbd'] < m['min_live_mbd']:
                    m['min_live_mbd'] = s['min_live_mbd']
            m['nevt'] += s['nevt']
            if m['ps_phot'] <= 0 and s['ps_phot'] > 0:
                m['ps_phot'] = s['ps_phot']
            if m['ps_mbd'] <= 0 and s['ps_mbd'] > 0:
                m['ps_mbd'] = s['ps_mbd']
    return merged


def parse_joey_allz(path):
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                rn = int(parts[0])
            except ValueError:
                continue
            try:
                b30c = float(parts[3])
                b30u = float(parts[4])
            except ValueError:
                continue
            if rn < 47000 or rn > 55000:
                continue
            out[rn] = (b30c, b30u)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-glob',
                    default='/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root')
    ap.add_argument('--joey-list',
                    default='/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/lumi/trig29_30_31_allz.list')
    ap.add_argument('--out-list',
                    default='/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/lumi/lumi_from_slimtree_allz.list')
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--limit', type=int, default=0,
                    help='Limit number of input files (for smoke tests).')
    args = ap.parse_args()

    files = sorted(glob.glob(args.data_glob))
    if args.limit > 0:
        files = files[:args.limit]
    print(f'Found {len(files)} slimtree files')
    if not files:
        sys.exit(1)

    stats_list = []
    with ProcessPoolExecutor(max_workers=args.nproc) as pool:
        futures = {pool.submit(per_file_scan, f): f for f in files}
        for i, fut in enumerate(as_completed(futures)):
            fpath = futures[fut]
            s = fut.result()
            stats_list.append(s)
            print(f'[{i+1}/{len(files)}] {os.path.basename(fpath)} -> {len(s)} runs')

    merged = merge_stats(stats_list)
    joey = parse_joey_allz(args.joey_list)

    rows = []
    L_uc_slim_tot = L_corr_slim_tot = 0.0
    L_uc_joey_tot = L_corr_joey_tot = 0.0
    n_with_joey = 0

    for rn in sorted(merged):
        s = merged[rn]
        if s['max_live_mbd'] <= 0 or s['ps_phot'] <= 0:
            continue
        L_uc = s['max_live_mbd'] / s['ps_phot'] / (SIGMA_MBD_MB * PB_PER_MB)
        jc = joey.get(rn)
        if jc is not None and jc[1] > 0:
            fcorr = jc[0] / jc[1]
        else:
            fcorr = 1.0
        L_corr = L_uc * fcorr
        L_uc_slim_tot += L_uc
        L_corr_slim_tot += L_corr
        if jc is not None:
            L_uc_joey_tot += jc[1]
            L_corr_joey_tot += jc[0]
            n_with_joey += 1
        rows.append((rn, L_uc, L_corr,
                     (jc[0] if jc else 0.0), (jc[1] if jc else 0.0),
                     s['ps_phot'], s['ps_mbd'],
                     s['max_live_mbd'], s['min_live_mbd'],
                     s['max_raw_mbd'], s['max_scaled_mbd'],
                     s['max_live_phot'], s['nevt']))

    with open(args.out_list, 'w') as f:
        f.write('# Columns: RN L_UC_slim L_Corr_slim L_Corr_joey L_UC_joey '
                'PS_phot PS_mbd max_live_mbd min_live_mbd max_raw_mbd '
                'max_scaled_mbd max_live_phot nevt_slim\n')
        f.write(f'# sigma_MBD = {SIGMA_MBD_MB} mb, MBD_BIT={MBD_BIT}, PHOT_BIT={PHOT_BIT}\n')
        for row in rows:
            f.write('{} {:.7f} {:.7f} {:.7f} {:.7f} {:g} {:g} {} {} {} {} {} {}\n'
                    .format(*row))
        f.write('#\n')
        f.write(f'# L_UC_slim total   : {L_uc_slim_tot:.5f} pb^-1\n')
        f.write(f'# L_Corr_slim total : {L_corr_slim_tot:.5f} pb^-1\n')
        f.write(f'# L_UC_joey total   : {L_uc_joey_tot:.5f} pb^-1\n')
        f.write(f'# L_Corr_joey total : {L_corr_joey_tot:.5f} pb^-1\n')
        f.write(f'# Runs present in slimtree: {len(rows)}\n')
        f.write(f'# Runs with Joey entry   : {n_with_joey}\n')

    print()
    print(f'Per-run rows written : {len(rows)}')
    print(f'Runs with Joey entry : {n_with_joey}')
    print(f'L_UC_slim total   : {L_uc_slim_tot:.5f} pb^-1')
    print(f'L_Corr_slim total : {L_corr_slim_tot:.5f} pb^-1')
    print(f'L_UC_joey total   : {L_uc_joey_tot:.5f} pb^-1')
    print(f'L_Corr_joey total : {L_corr_joey_tot:.5f} pb^-1')
    if L_uc_joey_tot > 0:
        print(f'Ratio UC   slim/joey : {L_uc_slim_tot/L_uc_joey_tot:.5f}')
    if L_corr_joey_tot > 0:
        print(f'Ratio Corr slim/joey : {L_corr_slim_tot/L_corr_joey_tot:.5f}')


if __name__ == '__main__':
    main()
