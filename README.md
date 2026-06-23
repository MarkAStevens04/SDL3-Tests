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

copy files: `scp "<local_dir>" markymoo@balam.scinet.utoronto.ca:<path_to_remote>`
    - NOTE: `<path_to_remote>` will start with /home/markymoo/
    - NOTE: to copy path from local, hold "option" key when right clicking.

## wandb
Sync with offline:
```
conda activate rgfn
wandb sync /scratch/markymoo/wandb/wandb/offline-run-*
```

## Chimera-X

- Find center of non-standard residue: `measure center :85C`


### Running PDB cleanup
cctbx         (protonate, adjust sidechains, optimize H-Bonds)     https://cctbx.github.io/
meeko         (Box specification, flexible residues, etc.)         https://meeko.readthedocs.io/en/develop/colab_examples.html
AutoDock Vina (actual docking)

1. `conda activate pdb_cleanup`
2. 


### Ideas:
Let's make a super easy to use and modular script for running inference. I'd love for this to be super parallelized. Let's have our script queue an `anchor_dock.py` job on the cluster. Ideally, this would take up 1/4 of a node (32 CPUs) to its maximum ability. Let's try to run 4 jobs in parallel. The input could be a directory to a list of sdf files. We'll 