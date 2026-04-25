#!/usr/bin/env python3
"""
Compute direct/total and frag/total ratio as a function of pT for the Werner
Vogelsang calculation, and print the JETPHOX(iso)/Werner(inclusive) ratio.

For isolated JETPHOX vs inclusive Werner, naive expectation:
    JETPHOX(iso) / Werner(incl) ~ direct_fraction + (1 - iso_rejection)*frag_fraction
where iso_rejection depends on cone+eps.  At RHIC 200 GeV, isolation typically
removes 60-80% of the fragmentation component.  So we expect:
    ratio ~ f_dir + 0.3 * f_frag   (if 70% iso rejection on frag)
If Werner were also isolated, ratio would be ~1 up to PDF/FF differences.
"""
import numpy as np
WERNER = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/sphenix_nlo"

def load(scale_key):
    return np.loadtxt(f"{WERNER}/photons_newphenix_sc{scale_key}.dat")

def main():
    arr = load("1")
    pt   = arr[:,0]
    direct = arr[:,1]
    frag = arr[:,2]
    tot  = arr[:,3]

    print(f"{'pT [GeV]':>8} {'direct':>12} {'frag':>12} {'total':>12} {'f_dir':>8} {'f_frag':>8}")
    print("-"*70)
    mask = (pt>=5.0) & (pt<=40.0)
    for i in np.where(mask)[0]:
        fdir = direct[i]/tot[i]
        ffr  = frag[i]/tot[i]
        print(f"{pt[i]:>8.2f} {direct[i]:>12.4e} {frag[i]:>12.4e} {tot[i]:>12.4e} {fdir:>8.3f} {ffr:>8.3f}")

    print()
    print("Expected JETPHOX(iso)/Werner(incl) if JETPHOX's iso rejects")
    print(" 80% of frag (ATLAS-like):     ratio ~ f_dir + 0.2 f_frag")
    print(" 60% of frag (looser iso):     ratio ~ f_dir + 0.4 f_frag")
    print(" 30% of frag (weak iso cone):  ratio ~ f_dir + 0.7 f_frag")
    print()
    print(f"{'pT':>5}  {'f_dir':>6}  {'f_frag':>6}  {'80% rej':>8}  {'60% rej':>8}  {'30% rej':>8}")
    for i in np.where(mask)[0]:
        fd = direct[i]/tot[i]; fr = frag[i]/tot[i]
        print(f"{pt[i]:>5.1f}  {fd:>6.3f}  {fr:>6.3f}  {fd+0.2*fr:>8.3f}  {fd+0.4*fr:>8.3f}  {fd+0.7*fr:>8.3f}")

if __name__ == "__main__":
    main()
