#!/bin/bash
# ClusterSizeStudy input for the DI cluster-size appendix figure, rerun with the
# current selection from config_showershape.yaml (2026-09-10).
set -e
source /sphenix/u/shuhang98/setup.sh
cd /sphenix/user/shuhangli/ppg12/efficiencytool
root -l -b -q "ClusterSizeStudy.C+(\"$1\",\"\",1.0,\"/sphenix/user/shuhangli/ppg12/efficiencytool/config_showershape.yaml\")"
