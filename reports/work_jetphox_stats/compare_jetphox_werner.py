#!/usr/bin/env python3
"""
Compare JETPHOX NLO and Werner Vogelsang NLO photon cross-section
per PPG12 pT bin (8-36 GeV).

JETPHOX  (h_truth_pT): pb/GeV integrated over |y|<0.7 -> divide by deta=1.4 to get
                      pb/(GeV*dy) at midrapidity.
Werner   (photons_newphenix_sc*.dat): columns [pt, direct, frag, total]
         where (direct+frag) == total is the invariant cross-section
         E d3sigma/dp3 [pb/GeV^2]. The plot macro converts:
             yield_plot = (direct+frag) * 2pi * pt  -> pb/(GeV*dy)
"""
import os
import sys
import numpy as np
import uproot

NLO_DIR   = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO/rootFiles"
WERNER    = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/sphenix_nlo"
DETA      = 1.4

# PPG12 reco pT bins (from plotcommon.h / analysis config)
PT_EDGES  = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]

def load_jetphox(scale):
    path = f"{NLO_DIR}/jetPHOX_{scale}.root"
    with uproot.open(path) as f:
        h = f["h_truth_pT"]
        # h_truth_pT is already pb/GeV (divided by bin width) integrated |y|<0.7
        edges = h.axis().edges()
        vals  = h.values()
    # map to arrays
    lo = edges[:-1]; hi = edges[1:]; cen = 0.5*(lo+hi)
    # divide by deta
    vals_per_dy = vals / DETA
    return cen, lo, hi, vals_per_dy

def load_werner(scale):
    # scale key: sc1, sc05, sc2
    f = f"{WERNER}/photons_newphenix_sc{scale}.dat"
    data = np.loadtxt(f)
    pt     = data[:,0]
    yield_ = data[:,1]
    frag   = data[:,2]
    total  = data[:,3]
    # after the plot macro: yield * 2pi * pt -> pb/(GeV*dy)
    werner = (yield_ + frag) * 2*np.pi * pt
    werner_total = total * 2*np.pi * pt
    return pt, werner, werner_total, yield_, frag

def werner_bin_avg(pt, yvals, xmin, xmax):
    """Linear interpolation integral / (xmax-xmin), matches plot macro."""
    # Trapezoidal over finely-sampled interior points
    # Mirror the TGraph Eval then Integral behavior: piecewise-linear.
    xs = np.linspace(xmin, xmax, 2001)
    ys = np.interp(xs, pt, yvals)
    integral = np.trapz(ys, xs)
    return integral / (xmax - xmin)

def werner_point_at(pt, yvals, x):
    return float(np.interp(x, pt, yvals))

def main():
    # --- JETPHOX all 3 scales
    jp = {}
    for sc in ("05","10","20"):
        cen, lo, hi, vals = load_jetphox(sc)
        jp[sc] = dict(cen=cen, lo=lo, hi=hi, val=vals)

    # --- Werner all 3 scales
    we = {}
    for key, label in [("05","05"), ("1","10"), ("2","20")]:
        pt, w, wtot, yld, frag = load_werner(key)
        we[label] = dict(pt=pt, val=w, total=wtot, direct=yld, frag=frag)

    # Print Werner file header sanity check
    print("="*90)
    print("Werner file layout (photons_newphenix_sc1.dat)")
    print("="*90)
    pt, w10, wtot10, yld10, frag10 = load_werner("1")
    # Note: print RAW columns, not the pre-multiplied ones
    raw = np.loadtxt(f"{WERNER}/photons_newphenix_sc1.dat")
    print(f"{'pt':>7} {'col2=direct':>12} {'col3=frag':>12} {'col4=total':>12} {'col2+col3':>12} {'match?':>8}")
    for i in [0,1,2,3,4,5,6,7,8,9,37,38]:
        d = raw[i,1]; fr = raw[i,2]; tot = raw[i,3]; s = d+fr
        ok = "Y" if abs(s-tot)/tot < 1e-2 else "N"
        print(f"{raw[i,0]:>7.2f} {d:>12.4e} {fr:>12.4e} {tot:>12.4e} {s:>12.4e} {ok:>8}")
    print(f"Werner pT range: {pt.min():.2f} - {pt.max():.2f} GeV, {len(pt)} rows")
    print()

    # --- Per-bin table
    print("="*90)
    print("JETPHOX vs Werner (central scale, mu = pT)")
    print("Units: pb/(GeV*dy) at |y|<0.7 midrapidity")
    print("="*90)
    header = f"{'bin [GeV]':<12} {'JP pb/GeV':>12} {'We_pt=ctr':>12} {'We_binAvg':>12} {'JP/We(ctr)':>12} {'JP/We(avg)':>12} {'dJP_scale%':>12} {'dWe_scale%':>12}"
    print(header)
    print("-"*len(header))

    rows = []
    for i in range(len(PT_EDGES)-1):
        xlow = PT_EDGES[i]; xhigh = PT_EDGES[i+1]; xcen = 0.5*(xlow+xhigh)

        # JETPHOX bin value at scale 10
        # find matching bin in jp[10]
        idx = np.argmin(np.abs(jp["10"]["cen"]-xcen))
        # sanity
        jp_lo = jp["10"]["lo"][idx]; jp_hi = jp["10"]["hi"][idx]
        if abs(jp_lo-xlow)>0.5 or abs(jp_hi-xhigh)>0.5:
            # Maybe JETPHOX uses truth bins; try matching by [xlow,xhigh] containing xcen
            idxs = np.where((jp["10"]["lo"]<=xcen)&(jp["10"]["hi"]>xcen))[0]
            if len(idxs)==1: idx = idxs[0]
        jpv10 = jp["10"]["val"][idx]
        jpv05 = jp["05"]["val"][idx]
        jpv20 = jp["20"]["val"][idx]

        # Werner at bin center and bin-averaged
        we10_ctr = werner_point_at(we["10"]["pt"], we["10"]["val"], xcen)
        we10_avg = werner_bin_avg(we["10"]["pt"], we["10"]["val"], xlow, xhigh)
        we05_avg = werner_bin_avg(we["05"]["pt"], we["05"]["val"], xlow, xhigh)
        we20_avg = werner_bin_avg(we["20"]["pt"], we["20"]["val"], xlow, xhigh)

        ratio_ctr = jpv10/we10_ctr if we10_ctr>0 else np.nan
        ratio_avg = jpv10/we10_avg if we10_avg>0 else np.nan

        jp_scale_span = (jpv05-jpv20)/jpv10*100 if jpv10>0 else np.nan
        we_scale_span = (we05_avg-we20_avg)/we10_avg*100 if we10_avg>0 else np.nan

        print(f"{xlow:>4.0f}-{xhigh:<7.0f}{jpv10:>12.4e} {we10_ctr:>12.4e} {we10_avg:>12.4e} {ratio_ctr:>12.3f} {ratio_avg:>12.3f} {jp_scale_span:>12.1f} {we_scale_span:>12.1f}")
        rows.append((xlow, xhigh, jpv10, we10_ctr, we10_avg, ratio_ctr, ratio_avg,
                     jpv05, jpv20, we05_avg, we20_avg))

    # Summary: ratio of central values, trend over pT
    print()
    print("="*90)
    print("Ratio trend: JETPHOX(10) / Werner(10) at bin center")
    print("="*90)
    r = [row[5] for row in rows]
    print(f"  low bin (8-10 GeV):    ratio = {r[0]:.3f}")
    print(f"  mid bin (14-16 GeV):   ratio = {r[3]:.3f}")
    print(f"  high bin (28-32 GeV):  ratio = {r[10]:.3f}")
    print(f"  highest (32-36 GeV):   ratio = {r[11]:.3f}")
    print()

    # Scale envelope overlap: does the JETPHOX band overlap Werner band?
    print("="*90)
    print("Scale-envelope overlap per bin")
    print(" JP band: [min(JP05,JP20), max(JP05,JP20)]   (JP05 is mu=0.5pT, larger)")
    print(" We band: [min(We05,We20), max(We05,We20)]")
    print("="*90)
    print(f"{'bin':<12}{'JP lo':>11}{'JP hi':>11}{'We lo':>11}{'We hi':>11}{'overlap?':>11}")
    for row in rows:
        xlow,xhigh,jp10,wctr,wavg,rctr,ravg,jp05,jp20,we05,we20 = row
        jplo, jphi = min(jp05,jp20), max(jp05,jp20)
        welo, wehi = min(we05,we20), max(we05,we20)
        ov = "YES" if (jplo<=wehi and welo<=jphi) else "NO"
        print(f"{xlow:.0f}-{xhigh:<7.0f}{jplo:>11.4e}{jphi:>11.4e}{welo:>11.4e}{wehi:>11.4e}{ov:>11}")

if __name__ == "__main__":
    main()
