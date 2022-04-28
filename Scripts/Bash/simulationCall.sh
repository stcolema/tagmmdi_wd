# Script to submit all of the simulation jobs to the HPC
#
# Stephen
# 19/04/2022

# test_size=(0.1 0.3 0.5 0.7 0.9);
test_size=(0.1 0.3);
adjust_for_TL=true;
generate_indices=0;

if [ ${generate_indices} -eq 1 ]
then 
  bash generateCVTestIndicesAllDatasets.sh;
fi

for i in "${test_size[@]}"
do
  echo ${i}
  echo "bash cvHPCcall.sh -a ${adjust_for_TL} -i ${i}";
  bash cvHPCcall.sh -a ${adjust_for_TL} -n ${i};
  echo "bash knnCall.sh -a ${adjust_for_TL} -i ${i}";
  bash knnCall.sh -a ${adjust_for_TL} -i ${i};
done

