setwd("/Users/dilerhaji/Desktop/willistoni-cg/Bernard-repeatmasked/blast-cds/")

##################################################
All CAF1 CDS sequences blasted against ONT assembly 
##################################################


blst <- read.table("equi-all-cds.out")
caf <- read.csv("caf1-features.csv", stringsAsFactors = FALSE)
dwill <- read.table("dwil_caf1_assembly_report.txt")
dwill_order <- read.csv("dwill-scaffold-order.csv")
dwill_order$scaf <- dwill[match(dwill_order[,2], dwill[,5]), 7]





#save.image("blst.RData")

head(caf, 30)
head(blst, 30)

blst$product_accession <- unlist(lapply(lapply(str_split(blst[,1], "_"), "[", 4:5), function(x){paste(x[1], x[2], sep = "_")}))
blst$genomic_accession <- unlist(lapply(lapply(str_split(unlist(lapply(str_split(blst[,1], "[|]"), "[", 2)), "_"), "[", 1:2), function(x){paste(x[1], x[2], sep = "_")}))
#save.image("blst.RData")

## E-value cuttoff
blst <- blst[blst[,11] < 0.0000000001, ]

## Percent identity cutoff 
blst <- blst[blst[,3] > 90, ]





###############################
# Matching scaffold accessions from CAF1 assembly to those in blast output
# Blast output is tblastn CAF1 CDS against assembly contigs
###############################

cafsub <- caf[caf$product_accession %in% blst$product_accession,]
cafsub_list <- list()
for(i in unique(cafsub$genomic_accession)){
	x <- cafsub[cafsub$genomic_accession %in% i, ]
	cafsub_list[[i]] <- x[order(x$start),]
}




###############################
# Considering only CAF1 scaffolds with more than 100 CDSs
# Considering only CAF1 scaffolds with more than 1000 blast hits 
###############################

cafsub_list100 <- cafsub_list[as.character(names(which(unlist(lapply(cafsub_list, function(x){ dim(x)[1]})) > 100)))]
str(cafsub_list100)
length(cafsub_list100)

blst_list <- list()
for(i in unique(blst$genomic_accession)){
	x <- blst[blst$genomic_accession %in% i, ]
	blst_list[[i]] <- x[order(x[,9]),]
}

blst_list1000 <- blst_list[as.character(names(which(unlist(lapply(blst_list, function(x){ dim(x)[1]})) > 1000)))]





###############################
# Ordering CDS based on (1) order on CAF1 references  
# and (2) order of blast position on CAF1 reference 
###############################

blst_caf1_order <- list()
for(i in names(blst_list1000)[names(blst_list1000) %in% names(cafsub_list100)]){
	x <- cafsub_list100[[i]]
	y <- blst_list1000[[i]]
	y$caf1_start <- x[match(y$product_accession, x$product_accession), "start"]
	y$blast_start_order <- order(y[,9])
	y$caf1_start_order <- order(y$caf1_start)
	blst_caf1_order[[i]] <- y
	}
	
ggplot(data = blst_caf1_order[[3]], aes(x = blast_start_order, y = caf1_start_order)) +
	geom_point()
	
	




###############################
# Parsing into muller elements and ordering/orienting 
###############################


blst_caf1_order2 <- do.call("rbind", blst_caf1_order)
head(blst_caf1_order2)
table(dwill_order$scaf %in% blst_caf1_order2[,14])

blst_caf1_order2$muller <- dwill_order[match(blst_caf1_order2[,14], dwill_order$scaf), 1]
blst_caf1_order2$scaff_order <- dwill_order[match(blst_caf1_order2[,14], dwill_order$scaf), 4]


blst_caf1_ori <- data.frame()
for(i in unique(blst_caf1_order2$genomic_accession)){
	if( dwill_order[dwill_order$scaf == i, 3] == "-" ) {
		x <- blst_caf1_order2[blst_caf1_order2$genomic_accession %in% i, ]
		x2 <- x[order(x$caf1_start_order), ]
		x2$caf1_start_order <- rev(x2$caf1_start_order)
	} else {
		x <- blst_caf1_order2[blst_caf1_order2$genomic_accession %in% i, ]
		x2 <- x[order(x$caf1_start_order), ]
	}
	blst_caf1_ori <- rbind(blst_caf1_ori, x2)
}


colnames(blst_caf1_ori)
ggplot(blst_caf1_ori, aes(x = caf1_start, y = V9)) +
	geom_point(size = 0.3) +
	facet_wrap(~muller+genomic_accession, scales = "free")