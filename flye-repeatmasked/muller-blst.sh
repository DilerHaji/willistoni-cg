#!/bin/bash
#SBATCH --job-name=wtd
#SBATCH --account=ac_flyminer
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --time=72:00:00

module load blast

tblastx -db mullerA-caf1 -query D.equinoxialis.assembly.medaka.fasta -outfmt 6 -out equi-mullA.out

