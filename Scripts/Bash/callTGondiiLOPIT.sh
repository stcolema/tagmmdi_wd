#!/bin/bash
#!
#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################


#! Example call:
# sbatch Scripts/SLURM/callMDITgondii.sh --job-name=T_gondii --time=35:59:59 --array=1-1;

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J mdiTGondii

#! Which project should be charged:
##SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -A CWALLACE-SL3-CPU

#! How many whole nodes should be allocated?
#SBATCH --nodes=1

#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1

#! How much wallclock time will be required?
#SBATCH --time=03:25:00

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem": 
#SBATCH -p skylake-himem

#! The number of jobs to run
#SBATCH --array=1-5

# Rscript specific parameters
R=65000;
thin=50;
data_file="/home/sdc56/rds/hpc-work/tagmmdi/T_gondii/Data/TGondiiMDI_K_125_input.rds";
script_dir='/home/sdc56/rds/hpc-work/tagmmdi/Scripts/R/';
save_dir='/home/sdc56/rds/hpc-work/tagmmdi/T_gondii/Output/';
script="callTAGMTGondii.R";

# Random seed used to define the fold
seed=${SLURM_ARRAY_TASK_ID};

# Number of chains to run
n_chains=1;

# Receive parameters from command line
while getopts r:t:d:k:i:s: option; do
    case "${option}" in
        r) R=${OPTARG};;
        t) thin=${OPTARG};;
        d) data_file=${OPTARG};;
        k) K=${OPTARG};;
        # s) seed=${OPTARG};;
        i) n_chains=${OPTARG};;
        s) save_dir=${OPTARG};;
    esac
done

#! sbatch directives end here (put any additional directives above this line)

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

#! Full path to application executable: 
application="Rscript ${script_dir}${script} --data_file ${data_file} --R ${R} --thin ${thin} --seed ${seed} --n_chains ${n_chains} --save_dir ${save_dir}";

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
