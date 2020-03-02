#!/bin/bash
#SBATCH --job-name=wtd
#SBATCH --account=co_rosalind
#SBATCH --partition=savio
#SBATCH --qos=savio_lowprio
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --time=72:00:00


#../../wtdbg2/wtdbg2 -t 20 -i ../ont/D.equinoxialis.guppyHAC.passReads.fastq.gz -fo D.equinoxialis -L 5000

#../../wtdbg2/wtdbg2 -t 20 -i ../ont/D.insularis.allReads.guppy324.fastq.gz -fo D.insularis -L 5000

#../../wtdbg2/wtdbg2 -t 20 -i ../ont/D.paulistorum.L06.allReads.guppy324.fastq.gz -fo D.paulistorum.L06 -L 5000

#../../wtdbg2/wtdbg2 -t 20 -i ../ont/D.paulistorum.L12.guppyHAC.passReads.fastq.gz -fo D.paulistorum.L12 -L 5000

module load samtools

#../../wtdbg2/wtpoa-cns -t 16 -i D.equinoxialis.ctg.lay -fo D.equinoxialis.raw.fa
#../../wtdbg2/wtpoa-cns -t 16 -i D.insularis.ctg.lay -fo D.insularis.raw.fa
#../../wtdbg2/wtpoa-cns -t 16 -i D.paulistorum.L06.ctg.lay -fo D.paulistorum.L06.raw.fa
#../../wtdbg2/wtpoa-cns -t 16 -i D.paulistorum.L12.ctg.lay -fo D.paulistorum.L12.raw.fa

#../../minimap2/minimap2 -t16 -ax map-pb -r2k D.equinoxialis.raw.fa ../ont/D.equinoxialis.guppyHAC.passReads.fastq.gz | samtools sort -@4 >D.equinoxialis.bam
samtools view -F0x900 D.equinoxialis.bam | ../../wtdbg2/wtpoa-cns -t 16 -i D.equinoxialis.raw.fa -fo D.equinoxialis.cns.fa

#../../minimap2/minimap2 -t16 -ax map-pb -r2k D.insularis.raw.fa ../ont/D.insularis.allReads.guppy324.fastq.gz | samtools sort -@4 >D.insularis.bam
samtools view -F0x900 D.insularis.bam | ../../wtdbg2/wtpoa-cns -t 16 -i D.insularis.raw.fa -fo D.insularis.cns.fa

#../../minimap2/minimap2 -t16 -ax map-pb -r2k D.paulistorum.L06.raw.fa ../ont/D.paulistorum.L06.allReads.guppy324.fastq.gz | samtools sort -@4 >D.paulistorum.L06.bam
samtools view -F0x900 D.paulistorum.L06.bam | ../../wtdbg2/wtpoa-cns -t 16 -i D.paulistorum.L06.raw.fa -fo D.paulistorum.L06.cns.fa

#../../minimap2/minimap2 -t16 -ax map-pb -r2k D.paulistorum.L12.raw.fa ../ont/D.paulistorum.L12.guppyHAC.passReads.fastq.gz | samtools sort -@4 >D.paulistorum.L12.bam
samtools view -F0x900 D.paulistorum.L12.bam | ../../wtdbg2/wtpoa-cns -t 16 -i D.paulistorum.L12.raw.fa -fo D.paulistorum.L12.cns.fa

../		


minimap2 -x map-ont layout.fasta escherichia_coli_map006_r7_3.fastq > m_1.paf
racon escherichia_coli_map006_r7_3.fastq m_1.paf layout.fasta > consensus_1.fasta
	