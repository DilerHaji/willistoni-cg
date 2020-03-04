library("ggplot2")
library("ggpubr")
library(stringr)
library(tidyr)
library(seqinr)
library(data.table)

setwd("/Users/dilerhaji/Desktop/willistoni-cg/Bernard-repeatmasked")



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

setwd("/Users/dilerhaji/Desktop/willistoni-cg/Bernard-repeatmasked/minimap-asm20")

library(stringr)
library(tidyr)
library(ggplot2)
library(pheatmap)

eqframe <- list()
eqframe2 <- list()
eqframe_plots <- list()
mappings_heatmap <- list()

for(i in system("ls *csv", intern = TRUE)){
	
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

colnames(eqframe2[[1]])
### "mullerA"  "mullerB"  "mullerC"  "mullerD"  "mullerEF"  ###
names(eqframe2)
	
for(j in c("equi.csv", "insu.csv", "paul.csv", "trop.csv", "will.csv")) {

	minimap_gg_plots <- list()
	for(i in names(eqframe2)) {

		for(j in 1:5) {
			if(is.null(rownames(eqframe2[[i]][which(!is.na(eqframe2[[i]][,j])),]))){
				m <- names(which(!is.na(eqframe2[[i]][,j])))
				ctg <- unlist(lapply(str_split(m, ":"), "[", 1))
			} else {
				m <- eqframe2[[i]][which(!is.na(eqframe2[[i]][,j])),]
				ctg <- unlist(lapply(str_split(rownames(m), ":"), "[", 1))

			}
   
		csv <- read.csv(i, head = FALSE, stringsAsFactors = FALSE)
		csv$mull <- unlist(lapply(str_split(unlist(lapply(str_split(csv[,6], "-"), "[", 2)), "[.]"), "[", 1))
		csv2 <- csv[as.character(csv[,1]) %in% as.character(ctg) & as.character(csv$mull) %in% as.character(colnames(eqframe2[[i]])[j]) & csv$V11 > 1000,]

		g <- ggplot(csv2, aes(V8, V3, col = V1, size = V11)) +
			geom_point() +
			scale_color_brewer(palette = "Paired")
   
		name <- paste(i, j, sep = "_")
		minimap_gg_plots[[name]] <- g
   
		}
		}


	for(i in 1:length(names(minimap_gg_plots))) {
		ggsave(paste(names(minimap_gg_plots)[i], ".png", sep = ""), minimap_gg_plots[[i]])
		}
		}
		




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
- Pairwise species comparisions (all vs all) 
##################################################

setwd("/Users/dilerhaji/Desktop/willistoni-cg/Bernard-repeatmasked/minimap-asm20")


map <- list()
file_names <- system("ls *csv", intern = TRUE)
file_names <- file_names[unlist(lapply(str_split(file_names, "-"), "[", 2)) != "caf1.csv"]


dat <- read.csv(file_names[2], head = FALSE, stringsAsFactors = FALSE)
ref <- read.csv(paste(unlist(lapply(str_split(j, "-"), "[", 1)), "caf1.csv", sep = "-"),  head = FALSE, stringsAsFactors = FALSE)
ref$mull <- unlist(lapply(str_split(unlist(lapply(str_split(ref$V6, "-"), "[", 2)),"[.]"), "[", 1))

ref[match(dat[,1], ref[,1]), "mullmax"]


plot(dat[,3], dat[,8], size = dat[,11])

ggplot(dat, aes(x = V3, y = V8, size = V11)) +
	geom_point()




for(j in file_names){
		
	eq <- read.csv(j, head = FALSE, stringsAsFactors = FALSE)
	
	ref <- read.csv(paste(unlist(lapply(str_split(j, "-"), "[", 1)), "caf1.csv", sep = "-"),  head = FALSE, stringsAsFactors = FALSE)
	ref$mull <- unlist(lapply(str_split(unlist(lapply(str_split(ref$V6, "-"), "[", 2)),"[.]"), "[", 1))
	
	ref2 <- read.csv(paste(unlist(lapply(str_split(unlist(lapply(str_split(j, "-"), "[", 2)), "[.]" ), "[", 1)), "caf1.csv", sep = "-"),  head = FALSE, stringsAsFactors = FALSE)
	ref2$mull <- unlist(lapply(str_split(unlist(lapply(str_split(ref2$V6, "-"), "[", 2)),"[.]"), "[", 1))

	
	### Muller element associated with the query contig 
	### Associating each unique assembly contig to a muller element by finding muller element to which most bases mapped (max)
		m1 <- aggregate(ref[,11], by = list(ref[,1], ref$mull), function(x){ mean(x) })
		m1max_ref <- c()
		for(i in unique(m1[,1])){ 
			x <- m1[m1[,1] == i, ]
			m1max_ref <- rbind(m1max_ref, c(i, x[x[,3] == max(x[,3]), 2]))
			}
		ref$mullmax <- m1max_ref[match(ref[,1], m1max_ref[,1]), 2]
		eq$match1 <- match(eq[,1], ref[,1])
		eq$muller1 <- ref[match(eq[,1], ref[,1]), "mullmax"]
	
	
	### Muller element associated with the subject contig 
	### Associating each unique assembly contig to a muller element by finding muller element to which most bases mapped (max)
		m2 <- aggregate(ref2[,11], by = list(ref2[,1], ref2$mull), function(x){ mean(x) })
		m1max_ref2 <- c()
		for(i in unique(m2[,1])){ 
			x <- m2[m2[,1] == i, ]
			m1max_ref2 <- rbind(m1max_ref2, c(i, x[x[,3] == max(x[,3]), 2]))
			}
		ref2$mullmax <- m1max_ref2[match(ref2[,1], m1max_ref2[,1]), 2]
		eq$match2 <- match(eq[,6], ref2[,1])
		eq$muller2 <-  ref[match(eq[,6], ref2[,1]), "mullmax"]
	
	
	eq <- eq[!is.na(eq$match1) & !is.na(eq$match2), ]
	#eq <- eq[eq$muller1 == eq$muller2, ]
	eq$ident <- eq[,10]/eq[,11]
	eq[,1] <- paste(eq[,1], paste(eq[,2], "bp", sep = ""), sep = ":  ")
	
	eq$comparison <- rep(j, dim(eq)[1])
	map[[j]] <- eq
	}
save.image("3March20.RData")

mapframe <- data.frame(data.table::rbindlist(map))

ggplot(mapframe[mapframe[,"muller1"] == "mullerEF",], aes(x = V3, y = V8, col = muller2)) + 
	geom_point(size = 0.2) +
	facet_wrap(~comparison, scales = "free")










#########################################
## Getting chains of colinear mappings 
#########################################

file_names <- system("ls *csv", intern = TRUE)
file_names <- file_names[unlist(lapply(str_split(file_names, "-"), "[", 2)) != "caf1.csv"]

map_matching <- list()

for(j in file_names){
 dat <- read.csv(j, head = FALSE, stringsAsFactors = FALSE)
 ref <- read.csv(paste(unlist(lapply(str_split(j, "-"), "[", 1)), "caf1.csv", sep = "-"),  head = FALSE, stringsAsFactors = FALSE)
 ref$mull <- unlist(lapply(str_split(unlist(lapply(str_split(ref$V6, "-"), "[", 2)),"[.]"), "[", 1))
 ref2 <- read.csv(paste(unlist(lapply(str_split(unlist(lapply(str_split(j, "-"), "[", 2)), "[.]" ), "[", 1)), "caf1.csv", sep = "-"),  head = FALSE, stringsAsFactors = FALSE)
 ref2$mull <- unlist(lapply(str_split(unlist(lapply(str_split(ref2$V6, "-"), "[", 2)),"[.]"), "[", 1))

 query_muller <- function(x){
	 mapping1 <- x[c(1,3)]
	 query <- ref[ref[,1] %in% as.character(mapping1[1]),]
	 dis1 <- abs(query[,3] - as.numeric(mapping1[2]))
	 dis1min <- query[which(dis1 == min(dis1)),]
	 if(dim(dis1min)[1] == 1) {
		 return(c(dis1min[, "mull"], dis1min[, 8]))
	 } else {
		 return("NA")
	 }
	 }

 subject_muller <- function(x){
	 mapping2 <- x[c(6,8)]
	 subject <- ref2[ref2[,1] %in% as.character(mapping2[1]),]
	 dis2 <- abs(subject[,3] - as.numeric(mapping2[2]))
	 dis2min <- subject[which(dis2 == min(dis2)),]
	 if(dim(dis2min)[1] == 1) {
		 return(c(dis2min[, "mull"], dis2min[, 8]))
	 } else {
		 return("NA")
	 }
	 } 

 query_muller_out <- apply(head(dat, 1000), 1, query_muller)
 subject_muller_out <- apply(head(dat, 1000), 1, subject_muller)

 dat <- cbind(dat, query_muller_out, subject_muller_out)

 map_matching[[j]] <- dat

}




dat$m1 <- ref[match(dat[,1], ref[,1]), "mull"]
dat$m1 <- ref[match(dat[,1], ref[,1]), "V3"]


dat$m2 <- ref2[match(dat[,6], ref2[,1]), "mull"]






head(dat)

dat <- dat[order(dat[,1], dat[,6], dat[,3]),]
head(dat)

mapbases <- aggregate(dat[,10]/dat[,11], by = list(dat[,1], dat[,6]), mean)
mapbases[order(mapbases[,1], mapbases[,2], mapbases[,3]),]

for(i in unique())