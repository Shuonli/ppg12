#!/usr/bin/env python3
"""Fill @FOO@ placeholders in converted_photon_study.tex with inclusive
efficiency numbers from rootFiles/efficiencies_summary.json.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE.parent / "converted_photon_study.tex"
JSON_IN = HERE / "rootFiles" / "efficiencies_summary.json"


def fmt_pct(x):
    return f"{100*x:.1f}\\%" if x is not None else "---"


def fmt_int(x):
    return f"{int(x):,}" if x is not None else "---"


def fmt_ratio(num, den):
    if den is None or den == 0 or num is None:
        return "---"
    return f"{num / den:.3f}"


def fmt_ratio_round(num, den):
    if den is None or den == 0 or num is None:
        return "0"
    return f"{num / den:.2f}"


def main():
    with open(JSON_IN) as f:
        s = json.load(f)
    unc = s["inclusive"]["unconv"]
    cnv = s["inclusive"]["conv"]

    subs = {
        "@NINCUNC@": fmt_int(unc["n_truth"]),
        "@NINCCNV@": fmt_int(cnv["n_truth"]),
        "@EPSRECOUNC@": fmt_pct(unc["eps_reco"]),
        "@EPSRECOCNV@": fmt_pct(cnv["eps_reco"]),
        "@EPSIDUNC@": fmt_pct(unc["eps_id"]),
        "@EPSIDCNV@": fmt_pct(cnv["eps_id"]),
        "@EPSISOUNC@": fmt_pct(unc["eps_iso"]),
        "@EPSISOCNV@": fmt_pct(cnv["eps_iso"]),
        "@EPSALLUNC@": fmt_pct(unc["eps_all"]),
        "@EPSALLCNV@": fmt_pct(cnv["eps_all"]),
        "@EPSIDCONDUNC@": fmt_pct(unc["eps_id_cond"]),
        "@EPSIDCONDCNV@": fmt_pct(cnv["eps_id_cond"]),
        "@EPSISOCONDUNC@": fmt_pct(unc["eps_iso_cond"]),
        "@EPSISOCONDCNV@": fmt_pct(cnv["eps_iso_cond"]),
        "@RATIORECO@": fmt_ratio(cnv["eps_reco"], unc["eps_reco"]),
        "@RATIOID@": fmt_ratio(cnv["eps_id"], unc["eps_id"]),
        "@RATIOISO@": fmt_ratio(cnv["eps_iso"], unc["eps_iso"]),
        "@RATIOALL@": fmt_ratio(cnv["eps_all"], unc["eps_all"]),
        "@RATIOIDCOND@": fmt_ratio(cnv["eps_id_cond"], unc["eps_id_cond"]),
        "@RATIOISOCOND@": fmt_ratio(cnv["eps_iso_cond"], unc["eps_iso_cond"]),
        "@RATIORECOROUND@": fmt_ratio_round(cnv["eps_reco"], unc["eps_reco"]),
        "@RATIOIDROUND@": fmt_ratio_round(cnv["eps_id"], unc["eps_id"]),
        "@RATIOISOCONDROUND@": fmt_ratio_round(cnv["eps_iso_cond"], unc["eps_iso_cond"]),
        "@RATIOALLROUND@": fmt_ratio_round(cnv["eps_all"], unc["eps_all"]),
    }

    src = TEX.read_text()
    missing = []
    for k, v in subs.items():
        if k in src:
            src = src.replace(k, v)
        else:
            missing.append(k)
    TEX.write_text(src)
    print(f"Filled {len(subs) - len(missing)} placeholders in {TEX.name}")
    if missing:
        print(f"  Missing from tex (no-op): {missing}")


if __name__ == "__main__":
    main()
