setwd("/Users/dilerhaji/Desktop/willistoni-cg/Bernard-repeatmasked/blast-cds/")

##################################################
All CAF1 CDS sequences blasted against ONT assembly 
##################################################


blst <- read.table("equi-all-cds.out")
caf <- read.csv("caf1-features.csv", stringsAsFactors = FALSE)
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



cafsub <- caf[caf$product_accession %in% blst$product_accession,]

head(cafsub)
dim(blst)

cafsub_list <- list()
for(i in unique(cafsub$genomic_accession)){
	x <- cafsub[cafsub$genomic_accession %in% i, ]
	cafsub_list[[i]] <- x[order(x$start),]
}

cafsub_list100 <- cafsub_list[as.character(names(which(unlist(lapply(cafsub_list, function(x){ dim(x)[1]})) > 100)))]
str(cafsub_list100)
length(cafsub_list100)



blst_list <- list()
for(i in unique(blst$genomic_accession)){
	x <- blst[blst$genomic_accession %in% i, ]
	blst_list[[i]] <- x[order(x[,9]),]
}

blst_list1000 <- blst_list[as.character(names(which(unlist(lapply(blst_list, function(x){ dim(x)[1]})) > 1000)))]


blst_caf1_order <- list()
for(i in names(blst_list1000)[names(blst_list1000) %in% names(cafsub_list100)]){
	x <- cafsub_list100[[i]]
	y <- blst_list1000[[i]]
	y$caf1_start <- x[match(y$product_accession, x$product_accession), "start"]
	y$blast_start_order <- order(y[,9])
	y$caf1_start_order <- order(y$caf1_start)
	blst_caf1_order[[i]] <- y
	}
	

length(blst_caf1_order[[1]]$blast_start_order)
length(unique(blst_caf1_order[[1]][,9]))


colnames(blst_caf1_order[[1]])
length(blst_caf1_order)

ggplot(data = blst_caf1_order[[13]], aes(x = blast_start_order, y = caf1_start_order)) +
	geom_point()