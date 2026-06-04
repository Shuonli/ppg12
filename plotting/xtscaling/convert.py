#!/usr/bin/env python3
"""
Convert (isolated) direct-photon cross sections from several experiments to the
common x_T-scaling variables of Bock HP2018 Fig.1 (left panel):

    x-axis: x_T   = 2 pT / sqrt(s)
    y-axis: Y     = (sqrt(s)/GeV)^n * E d^3sigma/dp^3   [ pb * GeV^-2 * c^3 ]   with n = 4.5

E d^3sigma/dp^3 (invariant cross section, at y=0) is recovered per dataset from its
reported observable:
    invariant_Ed3sigma_dp3        : I = C
    dsigma_dET_eta_integrated     : I = (C / delta_eta) / (2*pi*pT)
    d2sigma_dET_deta_per_eta      : I = C / (2*pi*pT)
Cross-section magnitude prefix (mb/nb/pb/fb) is taken from the units string.
Relative uncertainties are preserved under the linear transform, so
    Y_unc = Y * (unc_abs / C).

Inputs:
  data/ppg12_raw.csv                  PPG12 nominal (this analysis)
  data/phenix200_fig10.csv            PHENIX inclusive-direct (PRD86 072008)
  data/phenix200_isoratio_fig13a.csv  PHENIX isolated/inclusive ratio
  data/external_datasets.json         workflow output (ATLAS/CMS/E706 ...), schema list (optional)

Outputs:
  converted/<id>.csv   columns: xT, Y, Y_lo, Y_hi   (Y_lo/hi = Y -/+ total unc = quad(stat,syst))
  converted/_overlap_table.txt   sanity dump for the universal-curve check
"""
import os, csv, json, math

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT  = os.path.join(HERE, "converted")
os.makedirs(OUT, exist_ok=True)
N_SCALE = 4.5
TWO_PI = 2.0 * math.pi

def xs_prefix_factor(units):
    u = units.lower()
    if "cm" in u: return 1e36     # 1 cm^2 = 1e36 pb  (CGS invariant cross sections, e.g. R110/NA24)
    if "mb" in u: return 1e9      # 1 mb = 1e9 pb
    if "nb" in u: return 1e3      # 1 nb = 1e3 pb
    if "fb" in u: return 1e-3     # 1 fb = 1e-3 pb
    if "pb" in u: return 1.0
    return None

def invariant_from_point(C, observable, delta_eta, pT):
    if observable == "invariant_Ed3sigma_dp3":
        return C
    if observable == "dsigma_dET_eta_integrated":
        return (C / delta_eta) / (TWO_PI * pT)
    if observable == "d2sigma_dET_deta_per_eta":
        return C / (TWO_PI * pT)
    raise ValueError("unknown observable " + str(observable))

def convert_dataset(ds):
    """ds: schema dict. returns list of (xT, Y, Ytot_lo, Ytot_hi, Ystat, Ysyst)."""
    S = float(ds["sqrt_s_GeV"])
    obs = ds["observable"]
    deta = float(ds.get("delta_eta", 1) or 1)
    pref = xs_prefix_factor(ds["units"])
    if pref is None:
        raise ValueError(f"{ds['dataset_id']}: cannot parse units '{ds['units']}'")
    utype = ds.get("uncertainty_type", "absolute")
    rows = []
    for p in ds["points"]:
        if float(p["central"]) <= 0:
            continue  # skip zero / upper-limit points (e.g. UA6 last bin)
        x = float(p["x_ctr"])
        C = float(p["central"]) * pref
        I = invariant_from_point(C, obs, deta, x)
        Y = (S ** N_SCALE) * I
        xT = 2.0 * x / S
        # uncertainties: support symmetric stat, (a)symmetric syst
        def absunc(v):
            if v is None: return 0.0
            v = float(v)
            return (abs(v) / 100.0 * float(p["central"]) * pref) if utype == "percent" else abs(v) * pref
        stat = absunc(p.get("stat"))
        # syst may be single 'syst' or lo/hi
        syst_lo = absunc(p.get("syst_lo", p.get("syst")))
        syst_hi = absunc(p.get("syst_hi", p.get("syst")))
        if "total" in p and p.get("stat") is None and p.get("syst") is None:
            tot = absunc(p.get("total")); stat = tot; syst_lo = syst_hi = 0.0
        # PHENIX convention (PRD86 072008): exclude points consistent with zero,
        # i.e. whose lower error bar reaches/crosses zero (central <= stat (+) syst_lo).
        if abs(C) <= math.hypot(stat, syst_lo):
            continue
        # convert abs uncertainty (in central-units*pref) to Y via relative scaling
        Cmag = abs(C) if C != 0 else 1.0
        Ystat = Y * (stat / Cmag)
        Ysyst_lo = Y * (syst_lo / Cmag)
        Ysyst_hi = Y * (syst_hi / Cmag)
        tot_lo = math.hypot(Ystat, Ysyst_lo)
        tot_hi = math.hypot(Ystat, Ysyst_hi)
        # clamp the lower bound to stay positive: points whose lower uncertainty exceeds
        # the central value (>100%, consistent with zero) must not plunge to the log-axis floor
        y_lo = max(Y - tot_lo, Y * 0.02)
        rows.append((xT, Y, y_lo, Y + tot_hi, Ystat, 0.5 * (Ysyst_lo + Ysyst_hi)))
    rows.sort(key=lambda r: r[0])
    return rows

# ---------- build PPG12 entry ----------
def load_ppg12():
    pts = []
    with open(os.path.join(DATA, "ppg12_raw.csv")) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            a = [x.strip() for x in line.split(",")]
            if a[0].lower().startswith("et"): continue
            lo, hi, ctr, raw, stat, slo, shi = map(float, a[:7])
            # restrict to the reported, trigger-plateau range E_T in [12,36] GeV:
            # drop the [8,10],[10,12] turn-on bins and the [36,45] overflow bin
            if lo < 12.0 or hi > 36.0:
                continue
            pts.append({"x_lo": lo, "x_hi": hi, "x_ctr": ctr, "central": raw,
                        "stat": stat, "syst_lo": slo, "syst_hi": shi})
    return {"dataset_id": "PPG12_200GeV", "experiment": "sPHENIX (this analysis)",
            "sqrt_s_GeV": 200.0, "isolated": True, "observable": "dsigma_dET_eta_integrated",
            "x_variable": "ET", "units": "pb/GeV", "eta_range": "|eta|<0.7", "delta_eta": 1.4,
            "uncertainty_type": "absolute", "points": pts}

# ---------- build PHENIX entries ----------
def load_csv_points(fn):
    rows = []
    with open(os.path.join(DATA, fn)) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            rows.append([x.strip() for x in line.split(",")])
    return rows

def load_phenix_incl():
    pts = []
    for a in load_csv_points("phenix200_fig10.csv"):
        ctr, val, stm, stp, sy = map(float, a[:5])
        pts.append({"x_ctr": ctr, "central": val, "stat": 0.5*(stm+stp), "syst": sy})
    return {"dataset_id": "PHENIX_200GeV", "experiment": "PHENIX (incl.-direct)",
            "sqrt_s_GeV": 200.0, "isolated": False, "observable": "invariant_Ed3sigma_dp3",
            "x_variable": "pT", "units": "pb/GeV^2", "eta_range": "|eta|<0.25", "delta_eta": 1.0,
            "uncertainty_type": "absolute", "points": pts}

def load_phenix_isolated():
    """PHENIX inclusive-direct x (isolated/inclusive ratio), matched by nearest pT."""
    ratio = []
    for a in load_csv_points("phenix200_isoratio_fig13a.csv"):
        pt, r, st, sp, sm = map(float, a[:5]); ratio.append((pt, r))
    incl = load_phenix_incl()
    pts = []
    for p in incl["points"]:
        pt = p["x_ctr"]
        rr = min(ratio, key=lambda z: abs(z[0]-pt))[1]
        pts.append({"x_ctr": pt, "central": p["central"]*rr,
                    "stat": p["stat"]*rr, "syst": p["syst"]*rr})
    d = dict(incl); d["dataset_id"] = "PHENIX_200GeV_isolated"
    d["experiment"] = "PHENIX (isolated)"; d["isolated"] = True; d["points"] = pts
    return d

def load_alice7tev():
    """ALICE 7 TeV isolated photon (arXiv:1906.01371). d2sigma/(dpT deta) [nb/GeV/c], |eta|<0.27."""
    pts = []
    for a in load_csv_points("alice7tev.csv"):
        ctr, val, stm, stp, sy = map(float, a[:5])
        pts.append({"x_ctr": ctr, "central": val, "stat": 0.5*(stm+stp), "syst": sy})
    return {"dataset_id": "ALICE_7TeV", "experiment": "ALICE (7 TeV)",
            "sqrt_s_GeV": 7000.0, "isolated": True, "observable": "d2sigma_dET_deta_per_eta",
            "x_variable": "pT", "units": "nb/GeV", "eta_range": "|eta|<0.27", "delta_eta": 0.54,
            "uncertainty_type": "absolute", "points": pts}

def load_alice13tev():
    """ALICE 13 TeV isolated photon (arXiv:2407.01165). d2sigma/(dpT deta) [nb/GeV/c], |eta|<0.67, charged-only iso."""
    pts = []
    for a in load_csv_points("alice13tev.csv"):
        ctr, val, stm, stp, sy = map(float, a[:5])
        pts.append({"x_ctr": ctr, "central": val, "stat": 0.5*(stm+stp), "syst": sy})
    return {"dataset_id": "ALICE_13TeV", "experiment": "ALICE (13 TeV)",
            "sqrt_s_GeV": 13000.0, "isolated": True, "observable": "d2sigma_dET_deta_per_eta",
            "x_variable": "pT", "units": "nb/GeV", "eta_range": "|eta|<0.67", "delta_eta": 1.34,
            "uncertainty_type": "absolute", "points": pts}

def load_phenix510_iso():
    """PHENIX 510 GeV isolated direct photon (PRL 130, 251901). HEPData ins2033856 Fig 1 iso column."""
    pts = []
    for a in load_csv_points("phenix510_fig1.csv"):
        ctr, val, stm, stp, sy = map(float, a[:5])
        pts.append({"x_ctr": ctr, "central": val, "stat": 0.5*(stm+stp), "syst": sy})
    return {"dataset_id": "PHENIX_510GeV_isolated", "experiment": "PHENIX iso. (510 GeV)",
            "sqrt_s_GeV": 510.0, "isolated": True, "observable": "invariant_Ed3sigma_dp3",
            "x_variable": "pT", "units": "pb/GeV^2", "eta_range": "|eta|<0.25", "delta_eta": 1.0,
            "uncertainty_type": "absolute", "points": pts}

def main():
    datasets = [load_ppg12(), load_phenix_incl(), load_phenix510_iso(),
                load_alice7tev(), load_alice13tev()]
    extfn = os.path.join(DATA, "external_datasets.json")
    if os.path.exists(extfn):
        ext = json.load(open(extfn))
        for d in ext:
            if d.get("status") == "unavailable" or not d.get("points"):
                print(f"  skip {d.get('dataset_id')}: status={d.get('status')}"); continue
            datasets.append(d)
    else:
        print("  (no external_datasets.json yet -- PHENIX + PPG12 only)")

    overlap = open(os.path.join(OUT, "_overlap_table.txt"), "w")
    for ds in datasets:
        try:
            rows = convert_dataset(ds)
        except Exception as e:
            print(f"  ERROR {ds['dataset_id']}: {e}"); continue
        fn = os.path.join(OUT, ds["dataset_id"] + ".csv")
        with open(fn, "w") as f:
            f.write("# %s  sqrt(s)=%g GeV  obs=%s  units=%s  deta=%g  iso=%s\n" %
                    (ds["dataset_id"], ds["sqrt_s_GeV"], ds["observable"], ds["units"],
                     ds.get("delta_eta", 1), ds.get("isolated")))
            f.write("# xT, Y=(sqrt_s)^4.5*Ed3sigma/dp3, Y_lo, Y_hi\n")
            for r in rows:
                f.write("%.6e,%.6e,%.6e,%.6e\n" % (r[0], r[1], r[2], r[3]))
        msg = "%-26s n=%2d  xT[%.4f..%.4f]  Y[%.3e..%.3e]" % (
            ds["dataset_id"], len(rows), rows[0][0], rows[-1][0], rows[-1][1], rows[0][1])
        print("  wrote", msg)
        overlap.write(msg + "\n")
        for r in rows:
            overlap.write("    xT=%.5f  Y=%.4e\n" % (r[0], r[1]))
    overlap.close()
    print("done ->", OUT)

if __name__ == "__main__":
    main()
