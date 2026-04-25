#!/bin/bash
# Histogram driver for the CT14nlo JETPHOX rerun.
# Produces rootFiles/jetPHOX_nlo_{05,10,20}.root from pawres/ggd(o)rhic_nlo_{05,10,20}.root.
# Run this AFTER hadd_nlo.sh has merged the per-segment condor outputs.

cd "$(dirname "$0")"

root -b -q 'MakeJetPHOXhisto.C("05","_nlo")' > logs/jetPHOX_nlo_05.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("10","_nlo")' > logs/jetPHOX_nlo_10.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("20","_nlo")' > logs/jetPHOX_nlo_20.log 2>&1 &

wait
echo "All done. Outputs:"
ls -lh rootFiles/jetPHOX_nlo_{05,10,20}.root 2>/dev/null
