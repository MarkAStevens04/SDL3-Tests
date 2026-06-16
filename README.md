# SDL3-Tests

# Useful Commands
## Balam login
- Command + shift + p -> `Remote SSH: Connect to host...`
- Type in Passkey from password manager
- Type `1`, approve duo push from phone

## Slurm
View my jobs: `squeue --me`
View when my jobs are expected to run `squeue --me --start`
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

Set user-level variables: `nano ~/.bashrc`

View my scratch files: `ls $SCRATCH`

## wandb
Sync with offline:
```
conda activate rgfn
wandb sync /scratch/markymoo/wandb/wandb/offline-run-*
```

## Chimera-X

- Find center of non-standard residue: `measure center :85C`