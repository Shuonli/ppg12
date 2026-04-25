#!/bin/bash
# Runner for NNPDF3.1 JETPHOX histograms.
# Produces rootFiles/jetPHOX_nnpdf_{05,10,20}.root from pawres/ggd(o)rhic_nnpdf_{05,10,20}.root.

cd "$(dirname "$0")"

root -b -q 'MakeJetPHOXhisto.C("05","_nnpdf")' > logs/jetPHOX_nnpdf_05.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("10","_nnpdf")' > logs/jetPHOX_nnpdf_10.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("20","_nnpdf")' > logs/jetPHOX_nnpdf_20.log 2>&1 &

wait
echo "All done. Outputs:"
ls -lh rootFiles/jetPHOX_nnpdf_{05,10,20}.root 2>/dev/null
