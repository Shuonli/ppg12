#!/usr/bin/env python3
"""Independent re-derivation of JETPHOX per-bin cross-section + stat uncertainty.

Reads h_truth_pT (TH1F) from the three jetPHOX scale-variation ROOT files,
extracts per-bin edges, content (d sigma/dpT in pb/GeV), Sumw2 errors, and
derives stat_rel and N_eff. Also computes scale envelope and fiducial integrals.
"""
import uproot
import numpy as np

FILES = {
    "05": "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_05.root",
    "10": "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_10.root",
    "20": "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_20.root",
}


def read_hist(path):
    """Return (low_edges, high_edges, widths, content[pb/GeV], error[pb/GeV])."""
    with uproot.open(path) as f:
        h = f["h_truth_pT"]
        edges = h.axis().edges()
        # uproot returns values (content) including under/overflow? Use values(flow=False)
        content = h.values(flow=False)
        # Sumw2-propagated stat error via errors()
        errors = h.errors(flow=False)
    low = edges[:-1]
    high = edges[1:]
    width = high - low
    return low, high, width, content, errors


def print_table(scale, low, high, width, content, errors):
    print(f"\n===== Scale mu = {scale} (jetPHOX_{scale}.root) =====")
    print(f"{'idx':>3}  {'low':>6} {'high':>6} {'width':>6}  "
          f"{'content[pb/GeV]':>16} {'error[pb/GeV]':>14}  "
          f"{'stat_rel[%]':>11}  {'N_eff':>12}")
    for i, (lo, hi, w, c, e) in enumerate(zip(low, high, width, content, errors), 1):
        if c > 0:
            sr = e / c * 100.0
            neff = (c / e) ** 2 if e > 0 else np.inf
        else:
            sr = float("nan")
            neff = 0.0
        print(f"{i:>3}  {lo:>6.2f} {hi:>6.2f} {w:>6.2f}  "
              f"{c:>16.6g} {e:>14.4g}  {sr:>11.3f}  {neff:>12.1f}")


def integral_width(low, high, width, content, errors, pt_min=8.0, pt_max=36.0):
    """Integral of d sigma/dpT * d pT over [pt_min, pt_max], using only bins fully inside."""
    total = 0.0
    var = 0.0
    for lo, hi, w, c, e in zip(low, high, width, content, errors):
        if lo >= pt_min - 1e-6 and hi <= pt_max + 1e-6:
            total += c * w
            var += (e * w) ** 2
    return total, np.sqrt(var)


def main():
    data = {}
    for scale, path in FILES.items():
        data[scale] = read_hist(path)

    # Print per-scale tables
    for scale in ("05", "10", "20"):
        low, high, width, content, errors = data[scale]
        print_table(scale, low, high, width, content, errors)

    # Envelope table (bin-by-bin)
    low, high, width, c10, e10 = data["10"]
    _, _, _, c05, _ = data["05"]
    _, _, _, c20, _ = data["20"]

    print("\n===== Envelope per bin (using mu=1 as nominal) =====")
    print(f"{'idx':>3}  {'low':>6} {'high':>6}  "
          f"{'c05':>10} {'c10':>10} {'c20':>10}  "
          f"{'envelope_abs':>12} {'envelope_rel[%]':>15}  "
          f"{'stat_rel[%]':>11}  {'env/stat':>9}")
    for i in range(len(low)):
        lo, hi = low[i], high[i]
        a, b, c = c05[i], c10[i], c20[i]
        env_abs = max(a, b, c) - min(a, b, c)
        env_rel = env_abs / b * 100.0 if b > 0 else float("nan")
        stat_rel = e10[i] / b * 100.0 if b > 0 else float("nan")
        ratio = env_rel / stat_rel if stat_rel > 0 else float("inf")
        print(f"{i+1:>3}  {lo:>6.2f} {hi:>6.2f}  "
              f"{a:>10.4g} {b:>10.4g} {c:>10.4g}  "
              f"{env_abs:>12.4g} {env_rel:>15.3f}  "
              f"{stat_rel:>11.3f}  {ratio:>9.2f}")

    # Is env/stat >= 3 everywhere in 8-36 GeV?
    print("\n===== env/stat >= 3 check in 8-36 GeV =====")
    bad = []
    for i in range(len(low)):
        lo, hi = low[i], high[i]
        if lo < 8.0 - 1e-6 or hi > 36.0 + 1e-6:
            continue
        b = c10[i]
        if b <= 0:
            continue
        env_abs = max(c05[i], c10[i], c20[i]) - min(c05[i], c10[i], c20[i])
        env_rel = env_abs / b * 100.0
        stat_rel = e10[i] / b * 100.0
        ratio = env_rel / stat_rel if stat_rel > 0 else float("inf")
        flag = "OK" if ratio >= 3.0 else "FAIL"
        print(f"  bin [{lo:.1f},{hi:.1f}]: env_rel={env_rel:.2f}%, "
              f"stat_rel={stat_rel:.3f}%, ratio={ratio:.2f} {flag}")
        if ratio < 3.0:
            bad.append((lo, hi, ratio))
    if bad:
        print(f"  -> FAIL: {len(bad)} bins have env/stat < 3")
    else:
        print("  -> OK: env/stat >= 3 everywhere in 8-36 GeV")

    # Integrals over 8-36 GeV
    print("\n===== Fiducial integral 8-36 GeV (pb) =====")
    ints = {}
    for scale in ("05", "10", "20"):
        lo, hi, w, c, e = data[scale]
        I, dI = integral_width(lo, hi, w, c, e)
        ints[scale] = I
        print(f"  sigma(mu={scale}) = {I:.4g} +/- {dI:.3g} pb")
    print(f"  ratio 05/10 = {ints['05']/ints['10']:.4f}")
    print(f"  ratio 20/10 = {ints['20']/ints['10']:.4f}")

    # Wave 1 claim verification
    print("\n===== Wave 1 claim verification (mu=1) =====")
    claims = [
        (8.0, 10.0, 1519.0, 0.14),
        (16.0, 18.0, 47.01, 0.82),
        (28.0, 32.0, 0.9396, 4.09),
        (32.0, 36.0, 0.3227, 6.80),
    ]
    low10, high10, w10, c10_arr, e10_arr = data["10"]
    for lo_c, hi_c, c_claim, sr_claim in claims:
        found = False
        for i in range(len(low10)):
            if abs(low10[i] - lo_c) < 1e-3 and abs(high10[i] - hi_c) < 1e-3:
                found = True
                c = c10_arr[i]
                sr = e10_arr[i] / c * 100.0 if c > 0 else float("nan")
                d_c_rel = abs(c - c_claim) / c_claim * 100.0
                d_sr_abs = abs(sr - sr_claim)
                tag_c = "MATCH" if d_c_rel < 0.5 else "MISMATCH"
                tag_sr = "MATCH" if d_sr_abs < 0.05 else "MISMATCH"
                print(f"  bin [{lo_c},{hi_c}]: content extr={c:.4g} claim={c_claim} "
                      f"(diff {d_c_rel:.2f}%) {tag_c}; "
                      f"stat_rel extr={sr:.3f}% claim={sr_claim}% "
                      f"(abs diff {d_sr_abs:.3f}pp) {tag_sr}")
                break
        if not found:
            print(f"  bin [{lo_c},{hi_c}]: NO SUCH BIN in file")


if __name__ == "__main__":
    main()
