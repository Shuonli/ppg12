#!/bin/bash
# Histogram driver for the CT18NLO JETPHOX rerun.
# Produces rootFiles/jetPHOX_ct18_{05,10,20}.root from pawres/ggd(o)rhic_ct18_{05,10,20}.root.
# Run this AFTER hadd_ct18.sh has merged the per-segment condor outputs.

cd "$(dirname "$0")"

root -b -q 'MakeJetPHOXhisto.C("05","_ct18")' > logs/jetPHOX_ct18_05.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("10","_ct18")' > logs/jetPHOX_ct18_10.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("20","_ct18")' > logs/jetPHOX_ct18_20.log 2>&1 &

wait
echo "All done. Outputs:"
ls -lh rootFiles/jetPHOX_ct18_{05,10,20}.root 2>/dev/null
