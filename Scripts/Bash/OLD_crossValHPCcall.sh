#!/usr/bin/bash
#
# Calls to the cross-validation of datasets from ``pRolocdata``.
#
# Parameters of cross-validation call
R=5000;
burn=1000;
thin=25;
test_size=0.75;
cat_threshold=10;

# Submit the jobs to the HPC
sbatch callTagmMDICrossVal.sh -d "'tan2009r1 tan2009r1goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold};
sbatch callTagmMDICrossVal.sh -d "'HEK293T2011 HEK293T2011goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold};
sbatch callTagmMDICrossVal.sh -d "'E14TG2aS1 E14TG2aS1goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold};
sbatch callTagmMDICrossVal.sh -d "'groen2014r1 groen2014r1goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold};
sbatch callTagmMDICrossVal.sh -d "'dunkley2006 dunkley2006goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold};


