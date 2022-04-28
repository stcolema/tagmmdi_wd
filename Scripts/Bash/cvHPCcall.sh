#!/usr/bin/bash
#
# Calls to the cross-validation of datasets from ``pRolocdata``.
# Example: bash cvHPCcall.sh -r 10000 -t 25 -b 1000 -k 100 -c 10 -n 0.75 -s 10 -i 5;
#
# Parameters of cross-validation call
R=15000;
burn=5000;
thin=50;
test_size=0.10;
cat_threshold=5;
n_chains=5;
n_seeds=10;
K=100;
adjust_for_TL=false;

# Receive parameters from command line
while getopts r:t:b:k:c:n:s:i:a: option; do
    case "${option}" in
        r) R=${OPTARG};;
        t) thin=${OPTARG};;
        b) burn=${OPTARG};;
        k) K=${OPTARG};;
        c) cat_threshold=${OPTARG};;
        n) test_size=${OPTARG};;
        s) n_seed=${OPTARG};;
        i) n_chains=${OPTARG};;
        a) adjust_for_TL=${OPTARG};;
    esac
done

# Submit the jobs to the HPC
# for ((ii=1; ii<=${n_seeds}; ii++));
#do
sbatch --job-name="tanCV" --array=1-${n_seeds} --time="03:45:00" callTagmMDICVSingleFold.sh -d "'tan2009r1 tan2009r1goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold} -k ${K} -i ${n_chains} -a ${adjust_for_TL};
sbatch --job-name="hekCV" --array=1-${n_seeds} --time="17:30:00" callTagmMDICVSingleFold.sh -d "'HEK293T2011 HEK293T2011goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold} -k ${K} -i ${n_chains} -a ${adjust_for_TL};
sbatch --job-name="e14CV" --array=1-${n_seeds} --time="14:30:00" callTagmMDICVSingleFold.sh -d "'E14TG2aS1 E14TG2aS1goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold} -k ${K} -i ${n_chains} -a ${adjust_for_TL};
sbatch --job-name="groenCV" --array=1-${n_seeds} --time="06:30:00" callTagmMDICVSingleFold.sh -d "'groen2014r1 groen2014r1goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold} -k ${K} -i ${n_chains} -a ${adjust_for_TL};
sbatch --job-name="dunkleyCV" --array=1-${n_seeds} --time="04:15:00" callTagmMDICVSingleFold.sh -d "'dunkley2006 dunkley2006goCC'" -r ${R} -b ${burn} -t ${thin} -n ${test_size} -c ${cat_threshold} -k ${K} -i ${n_chains} -a ${adjust_for_TL};
# done;

