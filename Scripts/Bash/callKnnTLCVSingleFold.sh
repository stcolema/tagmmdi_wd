#!/bin/bash
#!
#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################


#! Example call:
# sbatch --job-name="knnHEK" --array=1-1 --time="35:59:59" callKnnTLCVSingleFold.sh -d "'HEK293T2011 HEK293T2011goCC'" -t "/home/sdc56/rds/hpc-work/tagmmdi/CV_output/HEK293T2011/HEK293T2011_HEK293T2011goCC_R_15000_seed_1_nChains_5_testSize_75.rds" -c 5 -n 5;

#! sbatch directives begin here ###############################
#! Name of the job:
###SBATCH -J tagmmdiCVSingleFold

#! Which project should be charged:
##SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -A CWALLACE-SL2-CPU

#! How many whole nodes should be allocated?
#SBATCH --nodes=1

#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1

#! How much wallclock time will be required?
###SBATCH --time=00:03:00

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem": 
#SBATCH -p skylake-himem

home_dir="/home/stephen/Documents/PhD/tagmmdi_wd";
datasets="'tan2009r1 tan2009r1goCC'";
script_dir='${home_dir}/Scripts/R/';
save_dir='${home_dir}/CV_output/';
test_index_dir="$save_dir";

# The threshold of entries in the categorical dataset for a feature to be kept
cat_threshold=5;

# Random seed used to define the fold
seed=${SLURM_ARRAY_TASK_ID};

# File holding information about test/train split
test_indices_file='tan2009r1_tan2009r1goCC_R_15000_seed_1_nChains_5_testSize_75.rds';
test_size=0.30;

# Number of weights considered for each class in the transfer learner
num_weights=5;

# Number of weight combinations used in optimisation
number_weights_sampled=10000;

# Were all classes represented in the training data
adjust_for_TL=false;

# Receive parameters from command line
while getopts d:c:t:n:s:i:a: option; do
    case "${option}" in
        d) datasets=${OPTARG};;
        c) cat_threshold=${OPTARG};;
        t) test_index_dir=${OPTARG};;
        n) num_weights=${OPTARG};;
        s) number_weights_sampled=${OPTARG};;
        i) test_size=${OPTARG};;
        a) adjust_for_TL=${OPTARG};;
    esac
done

#! sbatch directives end here (put any additional directives above this line)

IFS='\ ' read -r main_data auxiliary_data <<< "$datasets";

adjust_for_TL=${adjust_for_TL^^};
echo "Adjusting test indices for transfer learner: ${adjust_for_TL}";

test_str=$(printf %.0f $(echo "${test_size} * 100" | bc -l) );

main_dataset="$( cut -d ' ' -f 1 <<< "$datasets" )";
aux_dataset="$( cut -d ' ' -f 2 <<< "$datasets" )";

main_dataset="${main_dataset:1}";
aux_dataset="${aux_dataset:0:-1}";

R=15000;
n_chains=5;

test_index_dir="${test_index_dir}test_${test_str}/";
save_dir="${save_dir}test_${test_str}/";
test_indices_file="${test_index_dir}${main_dataset}/${main_dataset}_seed_${seed}_testSize_${test_str}_trainingAdjustedForTL_${adjust_for_TL}.rds";

# echo "${test_str}";
# echo "${test_index_dir}${main_dataset}/";
# exit;

# echo $main_dataset;
# echo $aux_dataset;
# echo $test_indices_file;
# exit;


#! Number of nodes and tasks per node allocated by SLURM (do not change):
numtasks=$SLURM_NTASKS
numnodes=$SLURM_JOB_NUM_NODES
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! ############################################################
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

# Load an R module
module load gcc-7.2.0-gcc-4.8.5-pqn7o2k
module load r-4.0.2-gcc-5.4.0-xyx46xb
module load netcdf/4.4.1


application="Rscript ${script_dir}knnTLCVSingleLoop.R --datasets ${datasets} --seed ${seed} --test_indices ${test_indices_file} --save_dir ${save_dir} --categorical_column_threshold ${cat_threshold} --number_weights ${num_weights} --number_weights_sampled ${number_weights_sampled}";

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.


cd $workdir

## If not submitting as a job to the cluster, don't do all these things

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
#np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
# CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

#! (OMP_NUM_THREADS threads will be created):
CMD="$application";


###############################################################
### You should not have to change anything below this line ####
###############################################################

# cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
echo "Time allocated: ${SBATCH_TIMELIMIT}";

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

#! (OMP_NUM_THREADS threads will be created):
# CMD="{ time $application $options; } 2> ${time_file} ";
CMD="$application $options";

eval $CMD 
