#!/usr/bin/env bash
# Run eta migration study for single and double interaction MC.
# Usage: bash run_eta_migration.sh [config_single] [config_double]
#   config_single  — YAML config for single interaction (default: config_bdt_nom.yaml)
#   config_double  — YAML config for double interaction (default: config_showershape.yaml)

source /sphenix/u/shuhang98/setup.sh

CONFIG_SINGLE=${1:-config_bdt_nom.yaml}
CONFIG_DOUBLE=${2:-config_showershape.yaml}

# Extract var_type from the single-interaction config
VAR_TYPE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIG_SINGLE}')); print(c['output']['var_type'])")
VAR_TYPE_DOUBLE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIG_DOUBLE}')); print(c['output']['var_type'])")

echo "============================================================"
echo "[run_eta_migration] CONFIG_SINGLE   : ${CONFIG_SINGLE}"
echo "[run_eta_migration] CONFIG_DOUBLE   : ${CONFIG_DOUBLE}"
echo "[run_eta_migration] VAR_TYPE        : ${VAR_TYPE}"
echo "[run_eta_migration] VAR_TYPE_DOUBLE : ${VAR_TYPE_DOUBLE}"
echo "============================================================"

# ==================================================================
# Single interaction: run EtaMigrationStudy for signal MC samples
# ==================================================================
echo ""
echo "[run_eta_migration] === Single interaction ==="
echo ""

root -l -b -q 'EtaMigrationStudy.C("'"${CONFIG_SINGLE}"'", "photon5")' &
root -l -b -q 'EtaMigrationStudy.C("'"${CONFIG_SINGLE}"'", "photon10")' &
root -l -b -q 'EtaMigrationStudy.C("'"${CONFIG_SINGLE}"'", "photon20")' &

wait
echo "[run_eta_migration] Single interaction jobs finished."

# Merge single-interaction results
hadd -f results/eta_migration_signal_${VAR_TYPE}.root \
  results/eta_migration_photon5_${VAR_TYPE}.root \
  results/eta_migration_photon10_${VAR_TYPE}.root \
  results/eta_migration_photon20_${VAR_TYPE}.root

echo "[run_eta_migration] Merged single → results/eta_migration_signal_${VAR_TYPE}.root"

# ==================================================================
# Double interaction: run EtaMigrationStudy with do_double=true
# ==================================================================
echo ""
echo "[run_eta_migration] === Double interaction ==="
echo ""

root -l -b -q 'EtaMigrationStudy.C("'"${CONFIG_DOUBLE}"'", "photon5",  true)' &
root -l -b -q 'EtaMigrationStudy.C("'"${CONFIG_DOUBLE}"'", "photon10", true)' &
root -l -b -q 'EtaMigrationStudy.C("'"${CONFIG_DOUBLE}"'", "photon20", true)' &

wait
echo "[run_eta_migration] Double interaction jobs finished."

# Merge double-interaction results
hadd -f results/eta_migration_double_signal_${VAR_TYPE_DOUBLE}.root \
  results/eta_migration_double_photon5_${VAR_TYPE_DOUBLE}.root \
  results/eta_migration_double_photon10_${VAR_TYPE_DOUBLE}.root \
  results/eta_migration_double_photon20_${VAR_TYPE_DOUBLE}.root

echo "[run_eta_migration] Merged double → results/eta_migration_double_signal_${VAR_TYPE_DOUBLE}.root"

# ==================================================================
# Summary
# ==================================================================
echo ""
echo "============================================================"
echo "[run_eta_migration] Output files:"
echo "  Single: results/eta_migration_signal_${VAR_TYPE}.root"
echo "  Double: results/eta_migration_double_signal_${VAR_TYPE_DOUBLE}.root"
echo "============================================================"
echo ""
echo "Next: plot with"
echo "  cd ../plotting && root -l -b -q 'plot_eta_migration.C(\"../efficiencytool/results/eta_migration_signal_${VAR_TYPE}.root\", \"../efficiencytool/results/eta_migration_double_signal_${VAR_TYPE_DOUBLE}.root\")'"
