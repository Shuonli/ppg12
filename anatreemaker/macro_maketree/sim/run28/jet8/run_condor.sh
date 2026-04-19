#!/usr/bin/bash

export macropath=$(pwd)

# Input list file
export listfile="${macropath}/test.list"
# Path to the output directory
export TargetDir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet8/condorout"
# Number of condor jobs
total_jobs=2000
# Skip lines (for the case of resubmission or continuation)
skip_lines=0

total_line=$(wc -l < "${listfile}")
lines_per_job=$(( (total_line + total_jobs - 1) / total_jobs ))

submitted=0
for ((i=0; i<total_jobs; i++)); do
    start_line=$(( i * lines_per_job + 1 + skip_lines))
    end_line=$(( start_line + lines_per_job - 1 + skip_lines))

    # Skip if start_line exceeds total lines
    if (( start_line > total_line )); then
        break
    fi

    mkdir -p "${TargetDir}/OutDir$i"
    WorkDir="${TargetDir}/OutDir$i"

    cp "$macropath/CondorRunSim_new.sh" "${WorkDir}/CondorRunTC$i.sh"
    cp "$macropath/Fun4All_run_sim.C" "${WorkDir}/"
    chmod +x "${WorkDir}/CondorRunTC$i.sh"

    sed -n "${start_line},${end_line}p" "${listfile}" > "${WorkDir}/test.list"

    submitted=$((submitted + 1))
done

# Build one shared submit file and submit once (replaces per-iteration condor_submit loop)
cat > "${TargetDir}/ff.sub" << EOF
+JobFlavour                   = "workday"
initialdir                    = ${TargetDir}/OutDir\$(Process)
transfer_input_files          = ${TargetDir}/OutDir\$(Process)/CondorRunTC\$(Process).sh,${TargetDir}/OutDir\$(Process)/Fun4All_run_sim.C,${TargetDir}/OutDir\$(Process)/test.list
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

Queue ${submitted}
EOF

condor_submit "${TargetDir}/ff.sub"
