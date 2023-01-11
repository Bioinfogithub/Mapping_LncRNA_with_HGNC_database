setwd("E:\\mtech project papers & codes 2022-23\\main data\\delncrna of GSE74697")
deg<- read.csv("deseq2_result.csv" , sep = "," , header = TRUE)
data <- read.csv("Processed_FPKM_data.csv", sep = "," , header = TRUE)
head(data)
dim(data)

## Remove columns with more than 50% NA
#dat[, which(colMeans(!is.na(dat)) > 0.5)]

## Remove rows with more than 50% NA
#dat[which(rowMeans(!is.na(dat)) > 0.5), ]

## Remove columns and rows with more than 50% NA
data1<- data[which(rowMeans(!is.na(data)) > 0.6), which(colMeans(!is.na(data)) > 0.6)]
dim(data1) 

#HGNC LncRNA list
lncrna <- read.delim("RNA_long_non-coding.txt")
head(lncrna)
head(deg)
#mapping gene with HGNC lncrna database . to check total number of lncrna in our deg file
lncrna.deg <- lncrna[lncrna$symbol %in% deg$Gene, ]
deg[lncrna.deg, ]
dim(lncrna.deg)
dim(deg)

library("org.Hs.eg.db") # remember to install it if you don't have it already
# symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
# 
# # Install the package if you have not installed by running this command: 
BiocManager::install("EnsDb.Hsapiens.v79")

#library(EnsDb.Hsapiens.v79)
library(biomaRt)
library("biomaRt")
 
library("EnsDb.Hsapiens.v79")
#converting gene symbol into ensemble id
geneIDs2<- ensembldb::select(EnsDb.Hsapiens.v79, keys= deg$Gene, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
head(geneIDs2)
dim(geneIDs2)
write.csv(geneIDs2, "DEG_Ensembl_ID.csv")

#mapping Enseble id in our deg with HGNC LncRNA ensembl_gene_id database . to find out thetotal number of LncRNA in our deg file.
lncrna.deg1 <- lncrna[lncrna$ensembl_gene_id %in% geneIDs2$GENEID, ]

dim(lncrna.deg1)
write.csv(lncrna.deg1, "Mapped_LncRNA.csv", sep = "," , quote = FALSE)

head(lncrna[,1:6])
#getting index of all lncrna in deg  
indx = which(deg$Gene %in% lncrna$symbol)
length(indx)
lncrna.deg = deg[indx,]#first col is index and , means give other value also in deg like- p value, logfold change, adj.p value
dim(lncrna.deg)
write.csv(lncrna.deg, "LncRNA_with_all_values.csv")

