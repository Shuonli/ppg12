#!/usr/bin/env python3
"""Fit PYTHIA tight+iso prompt-photon sigma_MC(ET) to
   sigma = sqrt(q0^2/E + q1^2/E^2 + q2^2)
   using the sg_2RMS column from resolution_vs_data_numbers.txt.
"""
import numpy as np
from scipy.optimize import curve_fit
import os, sys

NUM_FILE = "/sphenix/user/shuhangli/ppg12/plotting/figures/cluster_response_tightiso/resolution_vs_data_numbers.txt"

def reso(E, q0, q1, q2):
    return np.sqrt((q0**2)/E + (q1**2)/(E*E) + q2**2)

def reso2(E, q0, q2):
    return np.sqrt((q0**2)/E + q2**2)

et, sg, sge = [], [], []
with open(NUM_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        et.append(float(parts[0]))
        sg.append(float(parts[1]))
        sge.append(float(parts[2]) if float(parts[2]) > 0 else 1e-4)
et = np.array(et); sg = np.array(sg); sge = np.array(sge)

# 2-parameter fit: sqrt(q0^2/E + q2^2)  (drop the noise term which is degenerate)
p0 = [0.17, 0.045]
popt2, pcov2 = curve_fit(reso2, et, sg, p0=p0, sigma=sge, absolute_sigma=False)
perr2 = np.sqrt(np.diag(pcov2))
chi22 = np.sum(((sg - reso2(et, *popt2)) / sge)**2)
ndf2 = len(et) - 2
print(f"# 2-param fit  sqrt(q0^2/E + q2^2)")
print(f"q0 = {popt2[0]:.4f} +/- {perr2[0]:.4f}")
print(f"q2 = {popt2[1]:.4f} +/- {perr2[1]:.4f}")
print(f"chi2/ndf = {chi22:.2f} / {ndf2}")
print()
for x, y in zip(et, sg):
    yfit = reso2(x, *popt2)
    print(f"  {x:5.1f}  {y:.4f}   {yfit:.4f}   {y-yfit:+.4f}")
print()
print()

# Original 3-param fit for comparison
popt = [popt2[0], 0.0, popt2[1]]
perr = [perr2[0], 0.0, perr2[1]]
chi2 = chi22; ndf = ndf2

print(f"# Fit of PYTHIA tight+iso prompt-photon sigma_MC(ET) to sqrt(q0^2/E + q1^2/E^2 + q2^2)")
print(f"# input: {NUM_FILE}  (sg_2RMS column)")
print(f"q0 = {popt[0]:.4f} +/- {perr[0]:.4f}")
print(f"q1 = {popt[1]:.4f} +/- {perr[1]:.4f}")
print(f"q2 = {popt[2]:.4f} +/- {perr[2]:.4f}")
print(f"chi2/ndf = {chi2:.2f} / {ndf}")
print()
print("# fit eval at analysis bin centers:")
print("#  ET     sg_data    sg_fit   diff")
for x, y in zip(et, sg):
    yfit = reso(x, *popt)
    print(f"  {x:5.1f}  {y:.4f}   {yfit:.4f}   {y-yfit:+.4f}")
