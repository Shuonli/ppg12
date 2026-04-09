#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Execute the ROOT command with the argument
root -b -q "syst_eres.C" &
root -b -q "syst_escale.C" &
root -b -q "syst_fit.C" &
root -b -q "syst_iso.C" &
root -b -q "syst_mbd.C" &
root -b -q "syst_nor.C" &
root -b -q "syst_ntf.C" &
root -b -q "syst_tight.C" &
root -b -q "syst_fudge.C" &

wait
root -b -q "syst_purity.C" &
root -b -q "syst_eff.C" &
wait
root -b -q "calcSyst.C"



