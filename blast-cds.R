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


cafsub <- caf[caf$product_accession %in% blst$product_accession,]
head(cafsub)

