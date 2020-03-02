library("ggplot2")
library("ggpubr")
library(stringr)
library(tidyr)
library(seqinr)

setwd("~Desktop/Bernard-repeatmasked")



##################################################
Histograms of percent identity per muller element
##################################################


##### mullA ####
  equi <- read.csv("equi-mullA.csv", head = FALSE)
  colnames(equi)
  equi$prop <- equi[,10]/equi[,11]
  equi_histA <- ggplot(equi, aes(x = prop)) + 
	  geom_histogram(bins = 100) +
	  theme_bw()
  ggsave("equiA_hist.png", equi_hist, scale = 0.8)
  eq <- aggregate(equi$prop, by = list(equi[,1]), mean)
  hist(eq[,2], breaks = 100)
  eqsub <- eq[eq[,2] >= 0.85, 1]
  write.table(eqsub, "eqsub.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##### mullB ####
  equi <- read.csv("equi-mullB.csv", head = FALSE)
  colnames(equi)
  equi$prop <- equi[,10]/equi[,11]
  equi_histB <- ggplot(equi, aes(x = prop)) + 
	  geom_histogram(bins = 100) +
	  theme_bw()
  ggsave("equiB_hist.png", equi_hist, scale = 0.8)
  eq <- aggregate(equi$prop, by = list(equi[,1]), mean)
  hist(eq[,2], breaks = 100)
  eqsub <- eq[eq[,2] >= 0.85, 1]
  write.table(eqsub, "eqsub.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##### mullC ####
   equi <- read.csv("equi-mullC.csv", head = FALSE)
   colnames(equi)
   equi$prop <- equi[,10]/equi[,11]
   equi_histC <- ggplot(equi, aes(x = prop)) + 
	   geom_histogram(bins = 100) +
	   theme_bw()
   ggsave("equiC_hist.png", equi_hist, scale = 0.8)
   eq <- aggregate(equi$prop, by = list(equi[,1]), mean)
   hist(eq[,2], breaks = 100)
   eqsub <- eq[eq[,2] >= 0.85, 1]
   write.table(eqsub, "eqsub.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##### mullD ####
   equi <- read.csv("equi-mullD.csv", head = FALSE)
   colnames(equi)
   equi$prop <- equi[,10]/equi[,11]
   equi_histD <- ggplot(equi, aes(x = prop)) + 
	   geom_histogram(bins = 100) +
	   theme_bw()
   ggsave("equiD_hist.png", equi_hist, scale = 0.8)
   eq <- aggregate(equi$prop, by = list(equi[,1]), mean)
   hist(eq[,2], breaks = 100)
   eqsub <- eq[eq[,2] >= 0.85, 1]
   write.table(eqsub, "eqsub.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##### mullEF ####
   equi <- read.csv("equi-mullEF.csv", head = FALSE)
   colnames(equi)
   equi$prop <- equi[,10]/equi[,11]
   equi_histEF <- ggplot(equi, aes(x = prop)) + 
	   geom_histogram(bins = 100) +
	   theme_bw()
   ggsave("equiEF_hist.png", equi_hist, scale = 0.8)
   eq <- aggregate(equi$prop, by = list(equi[,1]), mean)
   hist(eq[,2], breaks = 100)
   eqsub <- eq[eq[,2] >= 0.85, 1]
   write.table(eqsub, "eqsub.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

ggg <- ggarrange(equi_histA, equi_histB, equi_histC, equi_histD , equi_histEF, labels = "AUTO")
ggsave("equi_hist.png", ggg)








##################################################
Heatmap of the number of minimap mappings to each muller element per assembly contig per species 
against CAF1 ordered/oriented muller elements 
##################################################

setwd("/Users/dilerhaji/Desktop/ont-genomes/Bernard-repeatmasked/minimap-asm20")


library(stringr)
library(tidyr)
library(ggplot2)
library(pheatmap)

eqframe <- list()
eqframe2 <- list()
eqframe_plots <- list()
mappings_heatmap <- list()

for(i in system("ls", intern = TRUE)){
	eq <- read.csv(i, head = FALSE, stringsAsFactors = FALSE)
	eq$mull <- unlist(lapply(str_split(unlist(lapply(str_split(eq$V6, "-"), "[", 2)),"[.]"), "[", 1))
	eq$ident <- eq[,10]/eq[,11]
	eq[,1] <- paste(eq[,1], paste(eq[,2], "bp", sep = ""), sep = ":  ")
	
	eq2 <- aggregate(eq[,11], by = list(eq[,1], eq$mull), function(x){length(x[x > 0])})
	len <- aggregate(eq[,2], by = list(eq[,1], eq$mull), max)
	ident <-  aggregate(eq$ident, by = list(eq[,1], eq$mull), function(x){mean(x)})
	mapsize <- aggregate(eq[,11], by = list(eq[,1], eq$mull), mean)

	eqframe[[i]] <- data.frame(contigs = eq2[,1], 
			muller = eq2[,2], 
			len = len[,3], 
			mapings = eq2[,3]*ident[,3], 
			ident = ident[,3],
			mapsize = mapsize[,3])
			
	eqframe_plots[[i]] <- ggplot(eqframe[[i]], aes(x = mapings, y = ident, col = len)) + 
								geom_point(size = 3) +
								facet_wrap(~muller)

	eqframe2[[i]] <- eqframe[[i]][eqframe[[i]]$mapings > 50, c(1,2,4)]
	eqframe2[[i]] <- spread(eqframe2[[i]], key = 2, value = 3)
	rownames(eqframe2[[i]]) <- eqframe2[[i]][,1]
	eqframe2[[i]] <- eqframe2[[i]][order(-as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(eqframe2[[i]][,1], ":"), "[", 2)), "bp"), "[", 1)))), ]
	eqframe2[[i]] <- as.matrix(eqframe2[[i]][, 2:6])
	eqframe2[[i]][is.na(eqframe2[[i]])] <- NA
	mappings_heatmap[[i]] <- pheatmap(eqframe2[[i]], na_col = "white", color = c("white", "grey","yellow", "orange", "red", "red4"), breaks = c(0, 10, 100, 1000, 2000, 3000, 4000), fontsize = 4, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 20)
}

ggsave("equi_mappings_heatmap.png", mappings_heatmap[["equi.csv"]] )
ggsave("trop_mappings_heatmap.png", mappings_heatmap[["trop.csv"]] )
ggsave("insu_mappings_heatmap.png", mappings_heatmap[["insu.csv"]] )
ggsave("paul_mappings_heatmap.png", mappings_heatmap[["paul.csv"]] )
ggsave("will_mappings_heatmap.png", mappings_heatmap[["will.csv"]] )









##################################################
- Extracting contigs sequences based on heatmap above 
- Ordering contig sequences based on majority mapping orientation of all mappings and CAF1 ordering 
- Quantifying overlap between assembly contigs 
##################################################

names(eqframe2)
i = "will.csv"
eqframe2[[i]]

m <- eqframe2[[i]][which(!is.na(eqframe2[[i]][,5])),]
ctg <- unlist(lapply(str_split(rownames(m), ":"), "[", 1))
csv <- read.csv(i, head = FALSE, stringsAsFactors = FALSE)
csv2 <- csv[as.character(csv[,1]) %in% as.character(ctg),]
csv2$order <- order(-csv2[,2])

ggplot(csv2, aes(V8, V2)) +
	geom_point()







mullerA_contigs <- unlist(lapply(str_split(rownames(eqframe2[[i]])[!is.na(eqframe2[[i]][,1] > 10)], ":"), "[", 1))
mullerB_contigs <- unlist(lapply(str_split(rownames(eqframe2)[!is.na(eqframe2[,2] > 10)], ":"), "[", 1))
mullerC_contigs <- unlist(lapply(str_split(rownames(eqframe2)[!is.na(eqframe2[,3] > 10)], ":"), "[", 1))
mullerD_contigs <- unlist(lapply(str_split(rownames(eqframe2)[!is.na(eqframe2[,4] > 10)], ":"), "[", 1))
mullerEF_contigs <- unlist(lapply(str_split(rownames(eqframe2)[!is.na(eqframe2[,5] > 10)], ":"), "[", 1))

write.table(mullerA_contigs, "mullerA_contigs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mullerB_contigs, "mullerB_contigs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mullerC_contigs, "mullerC_contigs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mullerD_contigs, "mullerD_contigs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mullerEF_contigs, "mullerEF_contigs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

system("seqkit grep -f mullerA_contigs.txt D.equinoxialis.assembly.sm.fasta > equi-mullA.fasta")
system("seqkit grep -f mullerB_contigs.txt D.equinoxialis.assembly.sm.fasta > equi-mullB.fasta")
system("seqkit grep -f mullerC_contigs.txt D.equinoxialis.assembly.sm.fasta > equi-mullC.fasta")
system("seqkit grep -f mullerD_contigs.txt D.equinoxialis.assembly.sm.fasta > equi-mullD.fasta")
system("seqkit grep -f mullerEF_contigs.txt D.equinoxialis.assembly.sm.fasta > equi-mullEF.fasta")










##################################################
- Candidate coding sequence set for blast 
##################################################

caf <- read.csv("caf1-features.csv", stringsAsFactors = FALSE)

#### Getting Dwill-caf1 contigs names for each muller element based on Schaeffer et al. 2006 GENETICS
caf1A <- unlist(lapply(str_split(system("cat mullerA-caf1.fasta | grep -o '>.*'", intern = TRUE), ">"), "[", 2))
caf1B <- unlist(lapply(str_split(system("cat mullerB-caf1.fasta | grep -o '>.*'", intern = TRUE), ">"), "[", 2))
caf1C <- unlist(lapply(str_split(system("cat mullerC-caf1.fasta | grep -o '>.*'", intern = TRUE), ">"), "[", 2))
caf1D <- unlist(lapply(str_split(system("cat mullerD-caf1.fasta | grep -o '>.*'", intern = TRUE), ">"), "[", 2))
caf1EF <- unlist(lapply(str_split(system("cat mullerEF-caf1.fasta | grep -o '>.*'", intern = TRUE), ">"), "[", 2))

cafall <- list(caf1A, caf1A, caf1C, caf1D, caf1EF)

#### Getting just the mapped contigs (above) and CDS features 
caf1A_fts <- caf[caf[,7] %in% caf1A & caf[,1] == "CDS", ]
caf1B_fts <- caf[caf[,7] %in% caf1B & caf[,1] == "CDS", ]
caf1C_fts <- caf[caf[,7] %in% caf1C & caf[,1] == "CDS", ]
caf1D_fts <- caf[caf[,7] %in% caf1D & caf[,1] == "CDS", ]
caf1EF_fts <- caf[caf[,7] %in% caf1EF & caf[,1] == "CDS", ]

cafall_fts <- list(caf1A_fts, caf1B_fts, caf1C_fts, caf1D_fts, caf1EF_fts)

#### Plotting distribution of CDS features and lengths across each contig 
ggplot(caf1EF_fts, aes(x = end, y = feature_interval_length)) +
	geom_point() +
	facet_wrap(~genomic_accession, scale = "free_x")


##### For each muller element (i) and contig (j), get 10 randomly picked CDS features and extract the coding sequence 
for(i in 1:5){
	for(j in 1:length(cafall[[i]])) {
		scaff <- cafall[[i]][j]
		fts <- cafall_fts[[i]][cafall_fts[[i]][,7] == scaff,]
		fts_sub <- fts[round(runif(10, 1, dim(fts)[1])),]
		gene_order <- fts_sub[order(fts_sub$start), 15]	
		write.table(gene_order, "gene_order.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
		system(paste("touch ", paste(scaff[1], "-", i, sep = ""), sep = ""))
		system(paste("grep -A 1 -f gene_order.txt caf1-cds2.fna | grep -v -- '^--$'", " > ", paste(scaff[1], "-", i, sep = ""), sep = ""))
			}	
		}







##################################################
- Blasting to assembly contigs and extracting position information 
##################################################

caf1_cds <- read.fasta("caf1-cds-translated.fasta")

str(caf1_cds)

attr(caf1_cds[[2]], "name")
length(caf1_cds[[2]])

caf1_cds_len <- lapply(caf1_cds, length)
caf1_cds_att <- lapply(caf1_cds, function(x){as.character(attr(x, "name"))})

genes <- unlist(lapply(str_split(unlist(lapply(str_split(unlist(caf1_cds_att), "GeneID:"), "[", 2)), "]"), "[", 1))
names(genes) <- 1:length(genes)
length(genes) - length(unique(genes))

for(i in unique(genes)){

	genes[genes %in% i]

}

hist(table(genes), breaks = 100)