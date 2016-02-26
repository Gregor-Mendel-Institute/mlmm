#set your working directory
setwd("/home/visegura/Doc-INRA/postdoc/Github_code_gwas/")

#load the tutorial data for carrying out mlmm analysis
genot <- read.table("data/genot.txt", sep = "\t", header = T)
genot_mat <- as.matrix(genot[, 2:ncol(genot)])
rownames(genot_mat) <- genot$Ind_id

phenot <- read.table("data/phenot.txt", sep = "\t", header = T)

map <- read.table("data/map.txt", sep = "\t", header = T)

genot_imp <- genot_mat
average <- colMeans(genot_imp, na.rm = T)

for (i in 1:ncol(genot_imp)){
  genot_imp[is.na(genot_imp[,i]), i] <- average[i]
}

stdev <- apply(genot_imp, 2, sd)
genot_stand <- sweep(sweep(genot_imp, 2, average, "-"), 2, stdev, "/")
K_mat <- (genot_stand %*% t(genot_stand)) / ncol(genot_stand)

#load the mlmm function as well as the emma package (if it does not install with your current R version, just download and source it, as recommended on the emma website).
source("mlmm.r")
source("emma.r")

#perform mlmm (10 steps), it can take few minutes...
mygwas <- mlmm(Y = phenot$Phenot1, X = genot_imp, K = K_mat, nbchunks = 2, maxsteps = 10, thresh = 1.2 * 10^-5)

#display and plot the results
source("plot_mlmm.r")
#mlmm stepwise table
mygwas$step_table
#EBIC plot
plot_step_table(mygwas, "extBIC")
#mbonf criterion plot
plot_step_table(mygwas, "maxpval") #user define threshold if defined is automatically drawn in red
#% variance plot
plot_step_RSS(mygwas)
#1st mlmm step plot
plot_fwd_GWAS(mygwas, step = 1, snp_info = map, pval_filt = 0.1)
#2nd mlmm step plot
plot_fwd_GWAS(mygwas, step = 2, snp_info = map, pval_filt = 0.1)
#3rd mlmm step plot
plot_fwd_GWAS(mygwas, step = 3, snp_info = map, pval_filt = 0.1)
#QQplot 7 steps
qqplot_fwd_GWAS(mygwas, nsteps = 7)

#optimal step according to ebic plot
plot_opt_GWAS(mygwas, opt = "extBIC", snp_info = map, pval_filt = 0.1)
qqplot_opt_GWAS(mygwas, opt = "extBIC")
#optimal step according to mbonf plot
plot_opt_GWAS(mygwas, opt = "mbonf", snp_info = map, pval_filt = 0.1)
qqplot_opt_GWAS(mygwas, opt = "mbonf")
#optimal step according to user defined threshold
plot_opt_GWAS(mygwas, opt = "thresh", snp_info = map, pval_filt = 0.1)#dotted line correspond to the user defined threshold
qqplot_opt_GWAS(mygwas, opt = "thresh")

#plot a region
plot_fwd_region(mygwas, step = 4, snp_info = map, pval_filt = 0.1, chrom = 2, pos1 = 17000000, pos2 = 19000000)
plot_opt_region(mygwas, opt = "thresh", snp_info = map, pval_filt = 0.1, chrom = 2, pos1 = 17000000, pos2 = 19000000)

#retrieving pvals
#step 1
head(mygwas$pval_step[[1]]$out)
#step 2
head(mygwas$pval_step[[2]]$out)
#opt extBIC
head(mygwas$opt_extBIC$out)
#including SNP effects
mygwas$opt_extBIC$coef

############
##including PCs to the model
PCs <- read.table("data/PCs.txt", sep = "\t", header = T)
PC_mat <- as.matrix(PCs[, 2:ncol(PCs)])
rownames(PC_mat) <- PCs$Ind_id

source("mlmm_cof.r")
mygwas_cof <- mlmm_cof(Y = phenot$Phenot1, X = genot_imp, cofs = PC_mat, K = K_mat, nbchunks = 2, maxsteps = 10, thresh = 10^-5)

#mlmm stepwise table
mygwas_cof$step_table
#EBIC plot
plot_step_table(mygwas_cof, "extBIC")
#mbonf criterion plot
plot_step_table(mygwas_cof, "maxpval")
#% variance plot
plot_step_RSS_cof(mygwas_cof)
#1st mlmm step plot
plot_fwd_GWAS(mygwas_cof, step = 1, snp_info = map, pval_filt = 0.1)
#2nd mlmm step plot
plot_fwd_GWAS(mygwas_cof, step = 2, snp_info = map, pval_filt = 0.1)
#3rd mlmm step plot
plot_fwd_GWAS(mygwas_cof, step = 3, snp_info = map, pval_filt = 0.1)
#QQplot 7 steps
qqplot_fwd_GWAS(mygwas_cof, nsteps = 7)

#optimal step according to ebic plot
plot_opt_GWAS(mygwas_cof, opt = "extBIC", snp_info = map, pval_filt = 0.1)
qqplot_opt_GWAS(mygwas_cof, opt = "extBIC")
#optimal step according to mbonf plot
plot_opt_GWAS(mygwas_cof, opt = "mbonf", snp_info = map, pval_filt = 0.1)
qqplot_opt_GWAS(mygwas_cof, opt = "mbonf")
#optimal step according to user defined threshold
plot_opt_GWAS(mygwas_cof, opt = "thresh", snp_info = map, pval_filt = 0.1)
qqplot_opt_GWAS(mygwas_cof, opt = "thresh")

#plot a region
plot_fwd_region(mygwas_cof, step = 4, snp_info = map, pval_filt = 0.1, chrom = 2, pos1 = 17000000, pos2 = 19000000)
plot_opt_region(mygwas_cof, opt = "thresh", snp_info = map, pval_filt = 0.1, chrom = 2, pos1 = 17000000, pos2 = 19000000)

#retrieving pvals
#step 1
head(mygwas_cof$pval_step[[1]]$out)
#step 2
head(mygwas_cof$pval_step[[2]]$out)
#opt extBIC
head(mygwas_cof$opt_extBIC$out)
#including SNP effects
mygwas_cof$opt_extBIC$coef

