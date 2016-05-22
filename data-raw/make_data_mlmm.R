library(devtools)

setwd("<path_to>/mlmm")

## load all input files
genot <- read.table("misc/genot.txt", sep = "\t", header = T)
genot_mat <- as.matrix(genot[, 2:ncol(genot)])
rownames(genot_mat) <- genot$Ind_id

phenot <- read.table("misc/phenot.txt", sep = "\t", header = T)

map <- read.table("misc/map.txt", sep = "\t", header = T)

PCs <- read.table("misc/PCs.txt", sep = "\t", header = T)
PC_mat <- as.matrix(PCs[, 2:ncol(PCs)])
rownames(PC_mat) <- PCs$Ind_id

## impute the missing genotypes and calculate the kinship matrix
genot_imp <- genot_mat
average <- colMeans(genot_imp, na.rm = T)
for (i in 1:ncol(genot_imp))
  genot_imp[is.na(genot_imp[,i]), i] <- average[i]
stdev <- apply(genot_imp, 2, sd)
genot_stand <- sweep(sweep(genot_imp, 2, average, "-"), 2, stdev, "/")
K_mat <- (genot_stand %*% t(genot_stand)) / ncol(genot_stand)

## format the data for the examples and save them
example_data <- list(X=genot_imp,
                     Y=phenot$Phenot1,
                     K=K_mat,
                     snp_info=map,
                     PC=PC_mat)
devtools::use_data(example_data, overwrite=TRUE)
