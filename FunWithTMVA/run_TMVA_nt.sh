#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if we have three arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <ivtx> <ieta> <iet>"
    exit 1
fi

# Run RecoEffCalculator for the various file types
root -l -b -q 'TMVA_train_nt.C('$1', '$2', '$3')'

