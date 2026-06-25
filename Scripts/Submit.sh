#!/bin/bash
sbatch <<EOT
#!/bin/sh
#SBATCH --job-name=S_$1                     # Job name 
#SBATCH --time=$2                           # Wall time
#SBATCH --mem-per-cpu=$3                    # Memory
#SBATCH --array 1-$4                        # nParams * nDistributions * nChunks - 1
#SBATCH --error=../Studies/S-${5}/Data/%a-${6}.ERR    # Errors output to approp. directory
#SBATCH --output=../Studies/S-${5}/Data/%a-${6}.OUT   # Simulation output to approp. directory

ml purge 
ml gcc/6.4.0 openmpi/3.0.0 boost/1.66.0

cd ../Studies/S-${5} && ./run.o null \${SLURM_ARRAY_TASK_ID} ${6}

# 1. Job name
# 2. Wall time
# 3. Memory limit
# 4. Array limits
# 5. Study name
# 6. Optional arguments list

# Example submit with optional arguments:

# for i in {0..10}; do
#   ./Submit MyJob 20:00:00 100MB 100 StudyName ${i}
# done


EOT
