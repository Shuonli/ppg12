#!/bin/bash

# Source the environment setup script
source /sphenix/u/shuhang98/setup.sh

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area: ${Cur_dir}"
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls
root '' > notes.log

# Run the ROOT macro with the file path as the second argument
root -b "Fun4All_runDST.C(\"inputdata.txt\",\"inputdatacalo.txt\",0)"

rm inputdata.txt
echo "DONE" >> DONE.txt

echo "JOB COMPLETE!"
