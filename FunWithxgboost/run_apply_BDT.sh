#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

INPUTFILE="$1"

WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
MACRO="$WORKDIR/apply_BDT.C"
CONFIG="config_nom.yaml"


cd "$WORKDIR"

echo "Running apply_BDT for data: $INPUTFILE"
root -l -b -q "$MACRO(\"$CONFIG\", \"data\", \"$INPUTFILE\")"
