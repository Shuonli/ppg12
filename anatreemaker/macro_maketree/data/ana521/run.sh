#!/bin/bash

export TargetDir="$PWD"/condorout


if compgen -G "${TargetDir}/OutDir*" > /dev/null; then
  echo "Removing existing OutDir* directories"
  rm -rf ${TargetDir}/OutDir*
else
  echo "No OutDir* directories found"
fi


i=0
echo "Creating a list of run numbers"
while read dir; do
  echo "Processing run number: $dir"
  li=$(printf "%04d" $i)

  rm inputdata.txt
  rm dst_jet*
  dstName="dstLists/dst_jet-000${dir}.list" 

  # creates a list of all files for a particular run
  if [ -f ${dstName} ]; then 
    cp ${dstName} inputdata.txt
  else
    echo "did not find ${dstName} getting from DB"
    CreateDstList.pl --tag ana521_2025p007_v001 --run ${dir}  DST_Jet
    cp dst_jet* inputdata.txt
    cp dst_jet* dstLists/.
  fi

  if [ ! -s inputdata.txt ]; then
    echo "inputdata.txt is empty, skipping to next iteration."
    continue
  fi
  #JETCALO-------------------------------------------------------------
  rm inputdatacalo.txt
  rm dst_jet*
  dstName="dstLists/dst_jetcalo-000${dir}.list"

  # creates a list of all files for a particular run
  if [ -f ${dstName} ]; then
    cp ${dstName} inputdatacalo.txt
  else
    echo "did not find ${dstName} getting from DB"
    CreateDstList.pl --tag ana521_2025p007_v001 --run ${dir}  DST_JETCALO
    cp dst_jet* inputdatacalo.txt
    cp dst_jet* dstLists/.
  fi

  if [ ! -s inputdatacalo.txt ]; then
    echo "inputdata.txt is empty, skipping to next iteration."
    continue
  fi

  j=10

  tot_files=$( cat inputdata.txt | wc -l )
  tot_calo_files=$( cat inputdatacalo.txt | wc -l )
  #check if the number of files in the two lists are the same
  if [ $tot_files -ne $tot_calo_files ]; then
    echo "The number of files in the two lists are not the same, skipping to next iteration."
    continue
  fi
  echo "total files: $tot_files"
  rem=$(( $tot_files%$j ))
  files_per_job=$(( $tot_files/$j ))
  njob=$j
  if [ $rem -ne 0 ]; then
    files_per_job=$(( $files_per_job+1 ))
  fi
  rem2=$(( $tot_files%$files_per_job ))
  njob=$(( $tot_files/$files_per_job ))
  if [ $rem2 -ne 0 ]; then
    njob=$(( ($tot_files/$files_per_job)+1 ))
  fi
  #echo "files per job: $files_per_job"
  #echo "njob: $njob"

  mkdir -p ${TargetDir}/OutDir$i
  export WorkDir="${TargetDir}/OutDir$i"
  echo "WorkDir:" ${WorkDir}

  cat > ${WorkDir}/ff.sub <<EOF
+JobFlavour                   = "workday"
Universe                      = vanilla
PeriodicHold                  = (NumJobStarts>=1 && JobStatus == 1)
Notification                  = Never
Priority                      = +80
request_memory                = 4GB
Notify_user                   = sl4859@columbia.edu

EOF

  pushd ${WorkDir}

  cp "$PWD"/../../Fun4All_runDST.C .
  cp "$PWD"/../../inputdata.txt .
  cp "$PWD"/../../inputdatacalo.txt .
  touch DONE.txt

  for((q=0;q<$njob;q++)); do
    start_file=$(( $q*$files_per_job+1 ))
    end_file=$(( $start_file+$files_per_job-1 ))
    #echo "start file: $start_file   end file: $end_file"

    sed -n $start_file\,${end_file}p inputdata.txt > ${WorkDir}/inputdata_$q.txt
    sed -n $start_file\,${end_file}p inputdatacalo.txt > ${WorkDir}/inputdatacalo_$q.txt

    cp "$PWD"/../../CondorRun.sh CondorRunJob${li}_$q.sh
    #cp "$PWD"/../../CondorRunWrap.sh CondorRunJob${li}_$q.sh
    sed -i "10 a cp ${WorkDir}/inputdata_$q.txt inputdata.txt" CondorRunJob${li}_$q.sh
    sed -i "10 a cp ${WorkDir}/inputdatacalo_$q.txt inputdatacalo.txt" CondorRunJob${li}_$q.sh
    sed -i "s/inputdata.txt/inputdata_$q.txt/g" CondorRunJob${li}_$q.sh
    sed -i "s/inputdatacalo.txt/inputdatacalo_$q.txt/g" CondorRunJob${li}_$q.sh
    echo "rm CondorRunJob${li}_$q.sh" >> CondorRunJob${li}_$q.sh
   #transfer_input_files          = ${WorkDir}/CondorRunJob${li}_$q.sh, ${WorkDir}/inputdata_$q.txt,${WorkDir}/Fun4All_EMCal.C

    cat >> ${WorkDir}/ff.sub <<EOF
Executable                    = ${WorkDir}/CondorRunJob${li}_$q.sh
Output                        =   ${WorkDir}/condor_${li}_$q.out
Error                         =  ${WorkDir}/condor_${li}_$q.errr
Log                           = /tmp/condor_sli_${li}_$q.log
Queue

EOF


    i=$((i+1))
  done

  echo "submitting some condor jobs"
  condor_submit ff.sub
  popd


done < runList.txt

# Set the directory where the files are located
file_directory="${TargetDir}/OutDir*/DONE.root"

while [ $(find condorout/OutDir* -name "DONE.txt" -print0 | xargs -0 cat | wc -l) -lt $((i)) ]; do
     current_file_count=$(find condorout/OutDir* -name "DONE.txt" -print0 | xargs -0 cat | wc -l)
    echo "Waiting for $((i)) files, currently $current_file_count"
    sleep 30  # Adjust the sleep duration as needed
done

