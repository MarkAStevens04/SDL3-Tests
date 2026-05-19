# SDL3-Tests

# Useful Commands
## Slurm
View my jobs: `sbatch --me`
View when my jobs are expected to run `sbatch --me --start`
See my utilization `sshare`
View queue `squeue`

Run a job:
```
cd $SCRATCH/rgfn_runs
sbatch /home/markymoo/projects/RGFN_Fork/RGFN-Fork/submit.sh
```
Note: You can set `#SBATCH --partition=debug` to get access to the debug node RIGHT NOW!

## Linux
Print current directory: `pwd`

## wandb
Sync with offline:
```
conda activate rgfn
wandb sync /scratch/markymoo/wandb/wandb/offline-run-*
```