#!/bin/bash
#SBATCH --job-name=mm
#SBATCH --account=ac_flyminer
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=2
#SBATCH --time=72:00:00


#../../../software/minimap2/minimap2 -cx asm10 all-caf1.fasta equi-hardmasked.fasta > equi.paf
#../../../software/minimap2/minimap2 -cx asm10 all-caf1.fasta insu-hardmasked.fasta > insu.paf
#../../../software/minimap2/minimap2 -cx asm10 all-caf1.fasta paul-hardmasked.fasta > paul.paf
#../../../software/minimap2/minimap2 -cx asm10 all-caf1.fasta trop-hardmasked.fasta > trop.paf
#../../../software/minimap2/minimap2 -cx asm10 all-caf1.fasta will-hardmasked.fasta > will.paf
../../../software/minimap2/minimap2 -cx asm10 equi-hardmasked.fasta equi-hardmasked.fasta > equi-equi.paf
../../../software/minimap2/minimap2 -cx asm10 insu-hardmasked.fasta equi-hardmasked.fasta > equi-insu.paf
../../../software/minimap2/minimap2 -cx asm10 paul-hardmasked.fasta  equi-hardmasked.fasta > equi-paul.paf
../../../software/minimap2/minimap2 -cx asm10 trop-hardmasked.fasta equi-hardmasked.fasta > equi-trop.paf
../../../software/minimap2/minimap2 -cx asm10 will-hardmasked.fasta equi-hardmasked.fasta > equi-will.paf
../../../software/minimap2/minimap2 -cx asm10 equi-hardmasked.fasta insu-hardmasked.fasta > insu-equi.paf
../../../software/minimap2/minimap2 -cx asm10 insu-hardmasked.fasta insu-hardmasked.fasta > insu-insu.paf
../../../software/minimap2/minimap2 -cx asm10 paul-hardmasked.fasta  insu-hardmasked.fasta > insu-paul.paf
../../../software/minimap2/minimap2 -cx asm10 trop-hardmasked.fasta insu-hardmasked.fasta > insu-trop.paf
../../../software/minimap2/minimap2 -cx asm10 will-hardmasked.fasta insu-hardmasked.fasta > insu-will.paf
../../../software/minimap2/minimap2 -cx asm10 equi-hardmasked.fasta paul-hardmasked.fasta > paul-equi.paf
../../../software/minimap2/minimap2 -cx asm10 insu-hardmasked.fasta paul-hardmasked.fasta > paul-insu.paf
../../../software/minimap2/minimap2 -cx asm10 paul-hardmasked.fasta  paul-hardmasked.fasta > paul-paul.paf
../../../software/minimap2/minimap2 -cx asm10 trop-hardmasked.fasta paul-hardmasked.fasta > paul-trop.paf
../../../software/minimap2/minimap2 -cx asm10 will-hardmasked.fasta paul-hardmasked.fasta > paul-will.paf
../../../software/minimap2/minimap2 -cx asm10 equi-hardmasked.fasta trop-hardmasked.fasta > trop-equi.paf
../../../software/minimap2/minimap2 -cx asm10 insu-hardmasked.fasta trop-hardmasked.fasta > trop-insu.paf
../../../software/minimap2/minimap2 -cx asm10 paul-hardmasked.fasta  trop-hardmasked.fasta > trop-paul.paf
../../../software/minimap2/minimap2 -cx asm10 trop-hardmasked.fasta trop-hardmasked.fasta > trop-trop.paf
../../../software/minimap2/minimap2 -cx asm10 will-hardmasked.fasta trop-hardmasked.fasta > trop-will.paf
../../../software/minimap2/minimap2 -cx asm10 equi-hardmasked.fasta will-hardmasked.fasta > will-equi.paf
../../../software/minimap2/minimap2 -cx asm10 insu-hardmasked.fasta will-hardmasked.fasta > will-insu.paf
../../../software/minimap2/minimap2 -cx asm10 paul-hardmasked.fasta  will-hardmasked.fasta > will-paul.paf
../../../software/minimap2/minimap2 -cx asm10 trop-hardmasked.fasta will-hardmasked.fasta > will-trop.paf
../../../software/minimap2/minimap2 -cx asm10 will-hardmasked.fasta will-hardmasked.fasta > will-will.paf
