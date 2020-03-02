#!/bin/bash
#SBATCH --job-name=wtdg2
#SBATCH --account=ac_flyminer
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=72:00:00

rclone copy Box:100x100/nanoporeReads/D.willistoni /global/scratch/diler/nanopore-genomes/bkim/ont/
rclone copy Box:100x100/nanoporeReads/D.equinoxialis /global/scratch/diler/nanopore-genomes/bkim/ont/
rclone copy Box:100x100/nanoporeReads/D.tropicalis /global/scratch/diler/nanopore-genomes/bkim/ont/
rclone copy Box:100x100/nanoporeReads/D.paulistorum.L06 /global/scratch/diler/nanopore-genomes/bkim/ont/
rclone copy Box:100x100/nanoporeReads/D.paulistorum.L12 /global/scratch/diler/nanopore-genomes/bkim/ont/
rclone copy Box:100x100/nanoporeReads/D.insularis /global/scratch/diler/nanopore-genomes/bkim/ont/
rclone copy Box:100x100/nanoporeReads/D.sturtevanti /global/scratch/diler/nanopore-genomes/bkim/ont/

rclone copy Box:100x100/illuminaForAssembly/D.willistoni /global/scratch/diler/nanopore-genomes/bkim/illumina/

rclone copy Box:100x100/assemblies/intermediate_steps/racon/D.equinoxialis* /global/scratch/diler/nanopore-genomes/bkim/flye-racon/
rclone copy Box:100x100/assemblies/intermediate_steps/racon/D.insularis* /global/scratch/diler/nanopore-genomes/bkim/flye-racon/
rclone copy Box:100x100/assemblies/intermediate_steps/racon/D.paulistorum* /global/scratch/diler/nanopore-genomes/bkim/flye-racon/
rclone copy Box:100x100/assemblies/intermediate_steps/racon/D.sturtevanti* /global/scratch/diler/nanopore-genomes/bkim/flye-racon/
rclone copy Box:100x100/assemblies/intermediate_steps/racon/D.tropicalis* /global/scratch/diler/nanopore-genomes/bkim/flye-racon/
rclone copy Box:100x100/assemblies/intermediate_steps/racon/D.willistoni* /global/scratch/diler/nanopore-genomes/bkim/flye-racon/

rclone copy Box:100x100/assemblies/intermediate_steps/madaka/D.equinoxialis* /global/scratch/diler/nanopore-genomes/bkim/flye-madaka/
rclone copy Box:100x100/assemblies/intermediate_steps/madaka/D.insularis* /global/scratch/diler/nanopore-genomes/bkim/flye-madaka/
rclone copy Box:100x100/assemblies/intermediate_steps/madaka/D.paulistorum* /global/scratch/diler/nanopore-genomes/bkim/flye-madaka/
rclone copy Box:100x100/assemblies/intermediate_steps/madaka/D.sturtevanti* /global/scratch/diler/nanopore-genomes/bkim/flye-madaka/
rclone copy Box:100x100/assemblies/intermediate_steps/madaka/D.tropicalis* /global/scratch/diler/nanopore-genomes/bkim/flye-madaka/
rclone copy Box:100x100/assemblies/intermediate_steps/madaka/D.willistoni* /global/scratch/diler/nanopore-genomes/bkim/flye-madaka/
