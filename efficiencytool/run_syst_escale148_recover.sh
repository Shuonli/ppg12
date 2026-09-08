#!/usr/bin/env bash
# run_syst_escale148_recover.sh  (2026-06-08)
#
# Race-free recovery for the escale148 syst run. The first Phase 2 (feeders +
# parents in one condor batch) hit the parent-vs-feeder merge race, breaking
# ~13 all-range merges and crashing CalcPhotonYield for energyscale148up
# (its Photon_final went missing). Feeders are all intact (213 keys).
#
#   Phase 2'  oneforall_syst_parents.sub  (21 bare parents, NO feeders -> no race)
#     -> wait
#   -> audit_merges.py
#   -> integrity re-check (print Photon_final integrals)
#   -> rebuild iso_generator variant
#   -> calc_syst_bdt.py

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u

LOG=logs/syst_escale148_recover.log
exec > >(tee -a "${LOG}") 2>&1

poll_cluster() {
    local cluster="$1"; local sleep_s="${2:-60}"
    sleep 10
    while condor_q "${cluster}" -nobatch 2>/dev/null | grep -qE "^${cluster}\."; do
        sleep "${sleep_s}"
    done
    set +e; condor_history -constraint "ClusterId==${cluster}" 2>/dev/null | head -3; set -e
}

echo "=========================================================="
echo "[recover] $(date -Is) START race-free parent re-run"
echo "[recover] parents: $(wc -l < syst_parents.list)"
echo "=========================================================="

P_OUT=$(condor_submit oneforall_syst_parents.sub)
echo "${P_OUT}"
P_CLUSTER=$(echo "${P_OUT}" | grep -oE "cluster [0-9]+" | grep -oE "[0-9]+" | tail -1)
echo "[recover] $(date -Is) parent cluster: ${P_CLUSTER}; polling every 60s..."
poll_cluster "${P_CLUSTER}" 60
echo "[recover] $(date -Is) parent batch cleared."

echo "[recover] $(date -Is) audit..."
set +e
python3 audit_merges.py > logs/audit_merges_escale148_recover.log 2>&1
echo "[recover] audit exit=$? (tail:)"; tail -3 logs/audit_merges_escale148_recover.log

echo "[recover] $(date -Is) integrity re-check (Photon_final integrals):"
python3 - <<'PYEOF'
import uproot, numpy as np
res="/sphenix/user/shuhangli/ppg12/efficiencytool/results/"
V=["nom","energyscale148up","energyscale148down","escale_nl","eres_smear_none","eres_smear_cE0p08",
   "noniso04","noniso10","npb03","npb07","purity_pade","purity_fit_ci_up","purity_fit_ci_down",
   "mc_purity_correction","mciso_no_shift","no_unfolding_reweighting","tightup_p05","ntdown_m10",
   "unfold_iter3","unfold_iter4","di_frac_fit"]
nom=None
for v in V:
    try:
        h=uproot.open(f"{res}Photon_final_bdt_{v}.root")["h_unfold_sub_result"].to_numpy()
        ig=float(np.nansum(h[0]))
        if v=="nom": nom=ig
        rel=f"{100*(ig/nom-1):+.1f}%" if nom else ""
        print(f"   {v:26s} integ={ig:8.4g}  {rel}")
    except Exception as e:
        print(f"   {v:26s} *** {e}")
PYEOF
set -e

echo "[recover] $(date -Is) rebuild iso_generator + aggregate..."
set +e
( cd ../plotting && root -l -b -q build_iso_generator_variant.C ) > logs/build_iso_generator_escale148_recover.log 2>&1
echo "[recover] iso_generator exit=$?"
( cd ../plotting && python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures --skip-missing ) > logs/calc_syst_bdt_escale148_recover.log 2>&1
echo "[recover] calc_syst_bdt exit=$?"
set -e

echo "=========================================================="
echo "[recover] $(date -Is) DONE (cluster ${P_CLUSTER})"
echo "=========================================================="
