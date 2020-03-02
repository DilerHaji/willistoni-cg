dwill <- read.table("Dwill_caf1_genome/ncbi-genomes-2020-01-27/GCF_000005925.1_dwil_caf1_assembly_report.txt")
dwill_order <- read.csv("Dwill_caf1_genome/dwill-scaffold-order.csv")


### Muller A

mullA <- dwill_order[dwill_order$muller == "A",]
write.table(dwill[dwill$V5 %in% as.character(mullA$gb), colnames(dwill) == "V7"], "mullerA", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("grep -A 1 -wFf mullerA Dwill_caf1_genome/GCF_000005925.1_dwil_caf1_genomic.fasta | sed '/^--$/d' > mullerA.fasta")
for(i in mullA$gb){
	temp <- as.character(dwill[dwill$V5 %in% as.character(i), colnames(dwill) == "V7"])
	write.table(temp, "temp", quote = FALSE, row.names = FALSE, col.names = FALSE)
	system("grep -A 1 -wFf temp | sed '/^--$/d' > temp.fasta")

}

### Muller B

mullB <- dwill_order[dwill_order$muller == "B",]
write.table(dwill[dwill$V5 %in% as.character(mullB$gb), colnames(dwill) == "V7"], "mullerB", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("grep -A 1 -wFf mullerB Dwill_caf1_genome/GCF_000005925.1_dwil_caf1_genomic.fasta | sed '/^--$/d' > mullerB.fasta")



### Muller C 

mullC <- dwill_order[dwill_order$muller == "C",]
write.table(dwill[dwill$V5 %in% as.character(mullC$gb), colnames(dwill) == "V7"], "mullerC", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("grep -A 1 -wFf mullerC Dwill_caf1_genome/GCF_000005925.1_dwil_caf1_genomic.fasta | sed '/^--$/d' > mullerC.fasta")



### Muller D 

mullD <- dwill_order[dwill_order$muller == "D",] 
write.table(dwill[dwill$V5 %in% as.character(mullD$gb), colnames(dwill) == "V7"], "mullerD", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("grep -A 1 -wFf mullerD Dwill_caf1_genome/GCF_000005925.1_dwil_caf1_genomic.fasta | sed '/^--$/d' > mullerD.fasta")



### Muller EF 

mullEF <- dwill_order[dwill_order$muller == "EF",] 
write.table(dwill[dwill$V5 %in% as.character(mullEF$gb), colnames(dwill) == "V7"], "mullerEF", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("grep -A 1 -wFf mullerEF Dwill_caf1_genome/GCF_000005925.1_dwil_caf1_genomic.fasta | sed '/^--$/d' > mullerEF.fasta")
