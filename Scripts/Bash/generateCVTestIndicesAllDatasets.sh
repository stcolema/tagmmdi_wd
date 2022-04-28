
datasets=('tan2009r1' 'HEK293T2011' 'groen2014r1' 'dunkley2006' 'E14TG2aS1')
test_size=(0.1 0.3 0.5 0.7 0.9);
adjust_for_TL=true;
n_seeds=10;

# Receive parameters from command line
while getopts a:n: option; do
    case "${option}" in
        a) adjust_for_TL=${OPTARG};;
        n) n_seeds=${OPTARG};;
    esac
done;
home_dir="/home/sdc56/rds/hpc-work/tagmmdi";
script_dir="${home_dir}/Scripts/R/";
gen_save_dir="${home_dir}/CV_output/test_";

for i in "${test_size[@]}"
do
  test_str=$(printf %.0f $(echo "${i} * 100" | bc -l) );
  save_dir="${gen_save_dir}${test_str}/";
  for dataset in "${datasets[@]}"
  do
    for ((j = 1; j < n_seeds + 1; j++ ));
    do
      Rscript ${script_dir}cvGenerateTestIndices.R --data ${dataset} --seed ${j} --test_size ${i} --save_dir ${save_dir} --adjust_for_TL ${adjust_for_TL} &
    done;
  done;
done;

