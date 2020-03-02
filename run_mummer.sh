#!/bin/bash
#SBATCH --job-name=mum
#SBatch --account=ac_flyminer
#SBATCH --account=ac_flyminer
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=72:00:00

## Command(s) to run:


MUMmer3/promer -p promer-equi-dmel-test Dmel.fasta D.equinoxialis.assembly.raconx4.fasta

MUMmer3/show-coords -r -c -l -L 100 -I 50 promer-equi-dmel.delta > promer-equi-dmel.coords


