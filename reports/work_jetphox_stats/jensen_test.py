#!/usr/bin/env python3
"""
Quantify the Werner bin-average bias from TGraph linear interpolation on
a steeply falling spectrum (Jensen's inequality).

TGraph->Integral(xmin,xmax)/(xmax-xmin) does piecewise-linear interpolation
between tabulated points. For a convex function like a steeply falling spectrum,
the linear chord OVER-estimates the true integral compared to the true
power-law curve. Quantify by fitting a local power law to the Werner points
and comparing the true integral to the linear-interpolation integral.
"""
import numpy as np
WERNER = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/sphenix_nlo"

def werner_values():
    arr = np.loadtxt(f"{WERNER}/photons_newphenix_sc1.dat")
    pt = arr[:,0]; y = (arr[:,1]+arr[:,2]) * 2*np.pi * pt
    return pt, y

def linear_integral(pt, y, xmin, xmax):
    xs = np.linspace(xmin, xmax, 2001)
    ys = np.interp(xs, pt, y)
    return np.trapz(ys, xs) / (xmax - xmin)

def powerlaw_integral(pt, y, xmin, xmax):
    """Use local power law: between adjacent points in log-log."""
    # identify the subrange of pt points inside [xmin,xmax]
    # subdivide the bin into sub-intervals matching Werner tabulation grid
    # power-law integral between two points (x1,y1) and (x2,y2):
    # y = A x^n, n = log(y2/y1)/log(x2/x1), A = y1 / x1^n
    # int = A/(n+1)*(x2^(n+1)-x1^(n+1)) if n != -1
    grid = np.concatenate(([xmin], pt[(pt>xmin)&(pt<xmax)], [xmax]))
    total = 0.0
    for i in range(len(grid)-1):
        x1 = grid[i]; x2 = grid[i+1]
        y1 = np.exp(np.interp(np.log(x1), np.log(pt), np.log(y)))
        y2 = np.exp(np.interp(np.log(x2), np.log(pt), np.log(y)))
        if x1 == x2 or y1<=0 or y2<=0:
            continue
        n = np.log(y2/y1)/np.log(x2/x1)
        if abs(n+1) < 1e-6:
            seg = y1*x1**n * np.log(x2/x1)  # degenerate case
        else:
            seg = (y1/x1**n) / (n+1) * (x2**(n+1) - x1**(n+1))
        total += seg
    return total / (xmax-xmin)

def main():
    pt, y = werner_values()
    print("Werner bin-average: LINEAR interp (plot macro) vs true POWER-LAW integral")
    print("Jensen bias = (linear - powerlaw) / powerlaw, positive = overestimate")
    print()
    print(f"{'bin [GeV]':<12}{'linear':>12}{'power-law':>12}{'bias %':>9}")
    edges = [8,10,12,14,16,18,20,22,24,26,28,32,36]
    for i in range(len(edges)-1):
        a = edges[i]; b = edges[i+1]
        lin = linear_integral(pt, y, a, b)
        pwl = powerlaw_integral(pt, y, a, b)
        bias = (lin - pwl)/pwl * 100
        print(f"{a:>4d}-{b:<6d}{lin:>12.4e}{pwl:>12.4e}{bias:>9.2f}")

if __name__ == "__main__":
    main()
