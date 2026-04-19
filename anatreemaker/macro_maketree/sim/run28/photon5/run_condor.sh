#!/usr/bin/bash

export macropath=$(pwd)

# Path to the dst files
export dstpath="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/photon5"
# Path to the output directory
export TargetDir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/photon5/condorout"
# Total events for this condor. (1 line = 1000 events)
total_line=10000
#total_line=10031
# Events per job (1 line = 1000 events)
lines_per_job=10
# Skip lines (for the case of resubmission or continuation)
skip_lines=0

total_jobs=$(( (total_line + lines_per_job - 1) / lines_per_job ))

for ((i=0; i<total_jobs; i++)); do
    start_line=$(( i * lines_per_job + 1 + skip_lines))
    end_line=$(( start_line + lines_per_job - 1 + skip_lines))

    mkdir -p "${TargetDir}/OutDir$i"
    WorkDir="${TargetDir}/OutDir$i"

    cp "$macropath/CondorRunSim_new.sh" "${WorkDir}/CondorRunTC$i.sh"
    cp "$macropath/Fun4All_run_sim.C" "${WorkDir}/"
    chmod +x "${WorkDir}/CondorRunTC$i.sh"

    sed -n "${start_line},${end_line}p" $dstpath/dst_calo_cluster.list > "${WorkDir}/dst_calo_cluster.list"
    sed -n "${start_line},${end_line}p" $dstpath/dst_truth_jet.list   > "${WorkDir}/dst_truth_jet.list"
    sed -n "${start_line},${end_line}p" $dstpath/g4hits.list          > "${WorkDir}/g4hits.list"
    sed -n "${start_line},${end_line}p" $dstpath/dst_mbd_epd.list     > "${WorkDir}/dst_mbd_epd.list"
done

# Build one shared submit file and submit once (replaces per-iteration condor_submit loop)
cat > "${TargetDir}/ff.sub" << EOF
+JobFlavour                   = "workday"
initialdir                    = ${TargetDir}/OutDir\$(Process)
transfer_input_files          = ${TargetDir}/OutDir\$(Process)/CondorRunTC\$(Process).sh,${TargetDir}/OutDir\$(Process)/Fun4All_run_sim.C
Executable                    = ${TargetDir}/OutDir\$(Process)/CondorRunTC\$(Process).sh
request_memory                = 4GB
PeriodicHold                  = (NumJobStarts>=2 && JobStatus == 1)
Universe                      = vanilla
Notification                  = Never
Priority                      = +80
Output                        = test.out
Error                         = test.err
Log                           = /tmp/sli_\$(Process).log
Notify_user                   = sl4859@columbia.edu

Queue ${total_jobs}
EOF

condor_submit "${TargetDir}/ff.sub"
