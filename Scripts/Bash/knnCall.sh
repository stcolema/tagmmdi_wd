

test_index_dir="/home/sdc56/rds/hpc-work/tagmmdi/CV_output/";
n_seeds=10;
c=5;
n=5;
s=20000;
test_size=0.50;
adjust_for_TL=false;

# Receive parameters from command line
while getopts c:t:n:s:i:a: option; do
    case "${option}" in
        c) c=${OPTARG};;
        t) test_index_dir=${OPTARG};;
        n) n=${OPTARG};;
        s) s=${OPTARG};;
        i) test_size=${OPTARG};;
        a) adjust_for_TL=${OPTARG};;
    esac
done

# test_str=$(echo "scale=0; ${test_size}*100" | bc);
test_str=$(printf %.0f $(echo "${test_size} * 100" | bc -l) );
# test_index_dir="${test_index_dir}${test_str}/";

sbatch --job-name="knnTan" --array=1-${n_seeds} --time="12:00:00" callKnnTLCVSingleFold.sh -d "'tan2009r1 tan2009r1goCC'" -t "${test_index_dir}" -c ${c} -n ${n} -s ${s} -i ${test_size} -a ${adjust_for_TL};
sbatch --job-name="knnDunkley" --array=1-${n_seeds} --time="12:00:00" callKnnTLCVSingleFold.sh -d "'dunkley2006 dunkley2006goCC'" -t "${test_index_dir}" -c ${c} -n ${n} -s ${s} -i ${test_size} -a ${adjust_for_TL};
sbatch --job-name="knnGroen" --array=1-${n_seeds} --time="12:00:00" callKnnTLCVSingleFold.sh -d "'groen2014r1 groen2014r1goCC'" -t "${test_index_dir}" -c ${c} -n ${n} -s ${s} -i ${test_size} -a ${adjust_for_TL};
sbatch --job-name="knnE14TG2" --array=1-${n_seeds} --time="34:00:00" callKnnTLCVSingleFold.sh -d "'E14TG2aS1 E14TG2aS1goCC'" -t "${test_index_dir}" -c ${c} -n ${n} -s ${s} -i ${test_size} -a ${adjust_for_TL};
sbatch --job-name="knnHEK" --array=1-${n_seeds} --time="34:00:00" callKnnTLCVSingleFold.sh -d "'HEK293T2011 HEK293T2011goCC'" -t "${test_index_dir}" -c ${c} -n ${n} -s ${s} -i ${test_size} -a ${adjust_for_TL};
