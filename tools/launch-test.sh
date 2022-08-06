#!/bin/bash

### ####################################################################################################
## Bash file to run a test. Function CompareFiles compares the generated 
## output files with the reference output files cpied in the reference directory
## #####################################################################################################

  CompareFile()
{
  ## Function to compare files. In case of OUT_EQCM, banner and appendix are removed from the comparison
  local genfile="$1"
  local ref="$2"

  if [  -f ${fileref} ]; then

    if [ $genfile = "OUT_EQCM" ] ;then
      nl=$(wc -l ${genfile} | awk '{ print $1 }')
      nappex=12
      nbanner=20
      nlow="$((nl - nappex))"
      ntop="$((nlow - nbanner))"
      (head -n $nlow $genfile | tail -n $ntop) &> out1.txt
      (head -n $nlow $fileref | tail -n $ntop) &> out2.txt
      diff out1.txt out2.txt  &> diff.log
      rm out1.txt out2.txt
    else
      diff $genfile $fileref &> diff.log
    fi

    if [ ! -s diff.log ]; then
      echo  "SUCCESS !!! ${genfile} passed the test"
    else
      echo  "FAILURE !!! ${genfile} has NOT passed the test"
    fi
    rm diff.log
  fi
}

# # Define the list of output files
   declare -a FileArray=("OUT_EQCM" "RESTART/RECORD_MODELS" \\
                         "ANALYSIS_EQCM/FILTERED_MASS"    "ANALYSIS_EQCM/SPEC_MASS"  "ANALYSIS_EQCM/RAW_MASS"\\
                         "ANALYSIS_EQCM/FILTERED_CURRENT" "ANALYSIS_EQCM/SPEC_CURRENT" "ANALYSIS_EQCM/RAW_CURRENT"\\
                         "ANALYSIS_EQCM/CHARACTERIZATION" "ANALYSIS_EQCM/MASS_CALIBRATION" "ANALYSIS_EQCM/ELECTRO_DEPOSITION" \\
                         "ANALYSIS_EQCM/INTERCALATION_OX" "ANALYSIS_EQCM/INTERCALATION_RED" "ANALYSIS_EQCM/MASSOGRAM"\\
                         "ANALYSIS_EQCM/SELECTED_SOLUTIONS" \\
                         "ATOMISTIC_MODELS/pristine/MODEL_SUMMARY" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/MODEL_SUMMARY" "ATOMISTIC_MODELS/1cycle-oxidation/model1/MODEL_SUMMARY" \\
                         "ATOMISTIC_MODELS/1cycle-reduction/MODEL_SUMMARY" "ATOMISTIC_MODELS/1cycle-reduction/model1/MODEL_SUMMARY" \\
                         "ATOMISTIC_MODELS/pristine/input.cp2k" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/input.cp2k" "ATOMISTIC_MODELS/1cycle-oxidation/model1/input.cp2k" \\
                         "ATOMISTIC_MODELS/1cycle-reduction/input.cp2k" "ATOMISTIC_MODELS/1cycle-reduction/model1/input.cp2k" \\
                         "ATOMISTIC_MODELS/pristine/INCAR" "ATOMISTIC_MODELS/pristine/KPOINTS" "ATOMISTIC_MODELS/pristine/POTCAR" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/INCAR"   "ATOMISTIC_MODELS/1cycle-oxidation/model1/INCAR" \\
                         "ATOMISTIC_MODELS/1cycle-reduction/INCAR"   "ATOMISTIC_MODELS/1cycle-reduction/model1/INCAR"  \\ 
                         "ATOMISTIC_MODELS/1cycle-oxidation/KPOINTS" "ATOMISTIC_MODELS/1cycle-oxidation/model1/KPOINTS" \\
                         "ATOMISTIC_MODELS/1cycle-reduction/KPOINTS" "ATOMISTIC_MODELS/1cycle-reduction/model1/KPOINTS" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/POTCAR"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/POTCAR"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/POTCAR"  "ATOMISTIC_MODELS/1cycle-reduction/model1/POTCAR" \\
                         "ATOMISTIC_MODELS/pristine/model.param" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/model.param"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/model.param"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/model.param"  "ATOMISTIC_MODELS/1cycle-reduction/model1/model.param" \\
                         "ATOMISTIC_MODELS/pristine/CI-castep.cell" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/CI-castep.cell"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/CI-castep.cell"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/CI-castep.cell"  "ATOMISTIC_MODELS/1cycle-reduction/model1/CI-castep.cell" \\
                         "ATOMISTIC_MODELS/pristine/CI-onetep.dat" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/CI-onetep.dat"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/CI-onetep.dat"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/CI-onetep.dat"  "ATOMISTIC_MODELS/1cycle-reduction/model1/CI-onetep.dat" \\
                         "ATOMISTIC_MODELS/pristine/hpc_script-vasp.sh" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/hpc_script-vasp.sh"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/hpc_script-vasp.sh"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/hpc_script-vasp.sh"  "ATOMISTIC_MODELS/1cycle-reduction/model1/hpc_script-vasp.sh" \\
                         "ATOMISTIC_MODELS/pristine/hpc_script-cp2k.sh" \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/hpc_script-cp2k.sh"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/hpc_script-cp2k.sh"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/hpc_script-cp2k.sh"  "ATOMISTIC_MODELS/1cycle-reduction/model1/hpc_script-cp2k.sh" \\
                         "ATOMISTIC_MODELS/pristine/hpc_script-castep.sh"\\
                         "ATOMISTIC_MODELS/1cycle-oxidation/hpc_script-castep.sh"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/hpc_script-castep.sh"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/hpc_script-castep.sh"  "ATOMISTIC_MODELS/1cycle-reduction/model1/hpc_script-castep.sh"  \\
                         "ATOMISTIC_MODELS/pristine/hpc_script-onetep.sh"  \\
                         "ATOMISTIC_MODELS/1cycle-oxidation/hpc_script-onetep.sh"  "ATOMISTIC_MODELS/1cycle-oxidation/model1/hpc_script-onetep.sh"  \\
                         "ATOMISTIC_MODELS/1cycle-reduction/hpc_script-onetep.sh"  "ATOMISTIC_MODELS/1cycle-reduction/model1/hpc_script-onetep.sh" )

  sep="/"	
  error="success"

  # Delete all output files before they are computed. This is because the test includes the
  # reference outputs files inside the testX.tar, which is also copied (and unpacked) to the "reference" directory.
  for i in "${FileArray[@]}"; do
    filename=${i}
    if [ -f ${filename} ]; then
      rm -rf $filename
    fi
  done

## Execute the job
$1

## Check if job has run ot not
if [ $? -ne 0 ]; then
  status="FAILED"
  echo "*** EXECUTION FAILED ***"  | tee -a diagnose.log 
  exit 125
else
  status="PASSED"
fi

## If the jobs for has executed, it is time to check the output with reference data
if [ "$status" = "PASSED" ]; then

  echo "*** Check for generated files against reference data ***"  | tee -a diagnose.log

  for i in "${FileArray[@]}"; do
    genfile=${i}
    fileref="$2$sep${i}"

    if [ -f ${genfile} ]; then
      CompareFile "$genfile" "$fileref" | tee -a diagnose.log
      check=$(awk '{w=$1} END{print w}' diagnose.log)
      if [ "$check" = "FAILURE" ] ; then
        error="fail"
      fi
    else
      if [ -f ${fileref} ]; then
        echo  "ERROR !! ${genfile} has NOT been generated"
        error="fail"
      fi
    fi
  done

  if [ $error = "fail" ] ; then
    echo "*** Test FAILED ****************************"  | tee -a diagnose.log
    exit 125
  else
    echo "*** Test was SUCCESSFUL ********************"  | tee -a diagnose.log
  fi  

fi
