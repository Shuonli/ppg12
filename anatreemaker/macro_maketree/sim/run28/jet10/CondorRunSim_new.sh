#!/usr/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

source /sphenix/u/shuhang98/setup.sh

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area:" ${Cur_dir}
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls

# Run the ROOT macro with the combined string as an argument
root  "Fun4All_run_sim.C"

echo "JOB COMPLETE!"
