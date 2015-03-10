##############################################################################################################################################
###MLMM - Multi-Locus Mixed Model
###SET OF FUNCTIONS TO CARRY GWAS CORRECTING FOR POPULATION STRUCTURE WHILE INCLUDING COFACTORS THROUGH A STEPWISE-REGRESSION APPROACH
#######
#
##note: require EMMA
#library(emma)
#source('emma.r')
#
##REQUIRED DATA & FORMAT
#
#PHENOTYPE - Y: a vector of length m, with names(Y)=individual names
#GENOTYPE - X: a n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names
#KINSHIP - K: a n by n matrix, with rownames(K)=colnames(K)=individual names
#each of these data being sorted in the same way, according to the individual name
#
##FOR PLOTING THE GWAS RESULTS
#SNP INFORMATION - snp_info: a data frame having at least 3 columns:
# - 1 named 'SNP', with SNP names (same as colnames(X)),
# - 1 named 'Chr', with the chromosome number to which belong each SNP
# - 1 named 'Pos', with the position of the SNP onto the chromosome it belongs to.
#######
#
##FUNCTIONS USE
#save this file somewhere on your computer and source it!
#source('path/mlmm.r')
#
###FORWARD + BACKWARD ANALYSES
#mygwas<-mlmm(Y,X,K,nbchunks,maxsteps)
#X,Y,K as described above
#nbchunks: an integer defining the number of chunks of X to run the analysis, allows to decrease the memory usage ==> minimum=2, increase it if you do not have enough memory
#maxsteps: maximum number of steps desired in the forward approach. The forward approach breaks automatically once the pseudo-heritability is close to 0,
#			however to avoid doing too many steps in case the pseudo-heritability does not reach a value close to 0, this parameter is also used.
#			It's value must be specified as an integer >= 3
#
###RESULTS
#
##STEPWISE TABLE
#mygwas$step_table
#
##PLOTS
#
##PLOTS FORM THE FORWARD TABLE
#plot_step_table(mygwas,type=c('h2','maxpval','BIC','extBIC'))
#
##RSS PLOT
#plot_step_RSS(mygwas)
#
##GWAS MANHATTAN PLOTS
#
#FORWARD STEPS
#plot_fwd_GWAS(mygwas,step,snp_info,pval_filt)
#step=the step to be plotted in the forward approach, where 1 is the EMMAX scan (no cofactor)
#snp_info as described above
#pval_filt=a p-value threshold for filtering the output, only p-vals below this threshold will be displayed in the plot
#
#OPTIMAL MODELS
#Automatic identification of the optimal models within the forwrad-backward models according to the extendedBIC or multiple-bonferonni criteria
#
#plot_opt_GWAS(mygwas,opt=c('extBIC','mbonf'),snp_info,pval_filt)
#snp_info as described above
#pval_filt=a p-value threshold for filtering the output, only p-vals below this threshold will be displayed in the plot
#
##GWAS MANHATTAN PLOT ZOOMED IN A REGION OF INTEREST
#plot_fwd_region(mygwas,step,snp_info,pval_filt,chrom,pos1,pos2)
#step=the step to be plotted in the forward approach, where 1 is the EMMAX scan (no cofactor)
#snp_info as described above
#pval_filt=a p-value threshold for filtering the output, only p-vals below this threshold will be displayed in the plot
#chrom is an integer specifying the chromosome on which the region of interest is
#pos1, pos2 are integers delimiting the region of interest in the same unit as Pos in snp_info
#
#plot_opt_region(mygwas,opt=c('extBIC','mbonf'),snp_info,pval_filt,chrom,pos1,pos2)
#snp_info as described above
#pval_filt=a p-value threshold for filtering the output, only p-vals below this threshold will be displayed in the plot
#chrom is an integer specifying the chromosome on which the region of interest is
#pos1, pos2 are integers delimiting the region of interest in the same unit as Pos in snp_info
#
##QQPLOTS of pvalues
#qqplot_fwd_GWAS(mygwas,nsteps)
#nsteps=maximum number of forward steps to be displayed
#
#qqplot_opt_GWAS(mygwas,opt=c('extBIC','mbonf'))
#
##############################################################################################################################################

mlmm_cof<-function(Y,X,cofs,K,nbchunks,maxsteps) {

n<-length(Y)
m<-ncol(X)

stopifnot(ncol(K) == n)
stopifnot(nrow(K) == n)
stopifnot(nrow(X) == n)
stopifnot(nrow(cofs) == n)
stopifnot(nbchunks >= 2)
stopifnot(maxsteps >= 3)

#INTERCEPT

Xo<-rep(1,n)

#K MATRIX NORMALISATION

K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
rm(K)

#step 0 : NULL MODEL

fix_cofs<-cbind(Xo,cofs)
rm(cofs)

addcof_fwd<-list()
addcof_fwd[[1]]<-'NA'

cof_fwd<-list()
cof_fwd[[1]]<-as.matrix(X[,colnames(X) %in% addcof_fwd[[1]]])

mod_fwd<-list()
mod_fwd[[1]]<-emma.REMLE(Y,cbind(fix_cofs,cof_fwd[[1]]),K_norm)

herit_fwd<-list()
herit_fwd[[1]]<-mod_fwd[[1]]$vg/(mod_fwd[[1]]$vg+mod_fwd[[1]]$ve)

RSSf<-list()
RSSf[[1]]<-'NA'

RSS_H0<-list()
RSS_H0[[1]]<-'NA'

df1<-1
df2<-list()
df2[[1]]<-'NA'

Ftest<-list()
Ftest[[1]]<-'NA'

pval<-list()
pval[[1]]<-'NA'

fwd_lm<-list()

cat('null model done! pseudo-h=',round(herit_fwd[[1]],3),'\n')

#step 1 : EMMAX

M<-solve(chol(mod_fwd[[1]]$vg*K_norm+mod_fwd[[1]]$ve*diag(n)))
Y_t<-crossprod(M,Y)
cof_fwd_t<-crossprod(M,cbind(fix_cofs,cof_fwd[[1]]))
fwd_lm[[1]]<-summary(lm(Y_t~0+cof_fwd_t))
Res_H0<-fwd_lm[[1]]$residuals
Q_<-qr.Q(qr(cof_fwd_t))

RSS<-list()
for (j in 1:(nbchunks-1)) {
X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% addcof_fwd[[1]]])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
rm(X_t)}
X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% addcof_fwd[[1]]])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof_fwd[[1]])))])
RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
rm(X_t,j)

RSSf[[2]]<-unlist(RSS)
RSS_H0[[2]]<-sum(Res_H0^2)
df2[[2]]<-n-df1-ncol(fix_cofs)-ncol(cof_fwd[[1]])
Ftest[[2]]<-(rep(RSS_H0[[2]],length(RSSf[[2]]))/RSSf[[2]]-1)*df2[[2]]/df1
pval[[2]]<-pf(Ftest[[2]],df1,df2[[2]],lower.tail=FALSE)
addcof_fwd[[2]]<-names(which(RSSf[[2]]==min(RSSf[[2]]))[1])
cof_fwd[[2]]<-cbind(cof_fwd[[1]],X[,colnames(X) %in% addcof_fwd[[2]]])
 colnames(cof_fwd[[2]])[ncol(cof_fwd[[2]])]<-addcof_fwd[[2]]
mod_fwd[[2]]<-emma.REMLE(Y,cbind(fix_cofs,cof_fwd[[2]]),K_norm)
herit_fwd[[2]]<-mod_fwd[[2]]$vg/(mod_fwd[[2]]$vg+mod_fwd[[2]]$ve)
rm(M,Y_t,cof_fwd_t,Res_H0,Q_,RSS)

cat('step 1 done! pseudo-h=',round(herit_fwd[[2]],3),'\n')

#FORWARD

for (i in 3:(maxsteps)) {
if (herit_fwd[[i-2]] < 0.01) break else {

M<-solve(chol(mod_fwd[[i-1]]$vg*K_norm+mod_fwd[[i-1]]$ve*diag(n)))
Y_t<-crossprod(M,Y)
cof_fwd_t<-crossprod(M,cbind(fix_cofs,cof_fwd[[i-1]]))
fwd_lm[[i-1]]<-summary(lm(Y_t~0+cof_fwd_t))
Res_H0<-fwd_lm[[i-1]]$residuals
Q_ <- qr.Q(qr(cof_fwd_t))

RSS<-list()
for (j in 1:(nbchunks-1)) {
X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[i-1]])])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
rm(X_t)}
X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[i-1]])])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof_fwd[[i-1]])))])
RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
rm(X_t,j)

RSSf[[i]]<-unlist(RSS)
RSS_H0[[i]]<-sum(Res_H0^2)
df2[[i]]<-n-df1-ncol(fix_cofs)-ncol(cof_fwd[[i-1]])
Ftest[[i]]<-(rep(RSS_H0[[i]],length(RSSf[[i]]))/RSSf[[i]]-1)*df2[[i]]/df1
pval[[i]]<-pf(Ftest[[i]],df1,df2[[i]],lower.tail=FALSE)
addcof_fwd[[i]]<-names(which(RSSf[[i]]==min(RSSf[[i]]))[1])
cof_fwd[[i]]<-cbind(cof_fwd[[i-1]],X[,colnames(X) %in% addcof_fwd[[i]]])
colnames(cof_fwd[[i]])[ncol(cof_fwd[[i]])]<-addcof_fwd[[i]]
mod_fwd[[i]]<-emma.REMLE(Y,cbind(fix_cofs,cof_fwd[[i]]),K_norm)
herit_fwd[[i]]<-mod_fwd[[i]]$vg/(mod_fwd[[i]]$vg+mod_fwd[[i]]$ve)
rm(M,Y_t,cof_fwd_t,Res_H0,Q_,RSS)}
cat('step ',i-1,' done! pseudo-h=',round(herit_fwd[[i]],3),'\n')}
rm(i)

##gls at last forward step
M<-solve(chol(mod_fwd[[length(mod_fwd)]]$vg*K_norm+mod_fwd[[length(mod_fwd)]]$ve*diag(n)))
Y_t<-crossprod(M,Y)
cof_fwd_t<-crossprod(M,cbind(fix_cofs,cof_fwd[[length(mod_fwd)]]))
fwd_lm[[length(mod_fwd)]]<-summary(lm(Y_t~0+cof_fwd_t))

Res_H0<-fwd_lm[[length(mod_fwd)]]$residuals
Q_ <- qr.Q(qr(cof_fwd_t))

RSS<-list()
for (j in 1:(nbchunks-1)) {
X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[length(mod_fwd)]])])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
rm(X_t)}
X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[length(mod_fwd)]])])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof_fwd[[length(mod_fwd)]])))])
RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
rm(X_t,j)

RSSf[[length(mod_fwd)+1]]<-unlist(RSS)
RSS_H0[[length(mod_fwd)+1]]<-sum(Res_H0^2)
df2[[length(mod_fwd)+1]]<-n-df1-ncol(fix_cofs)-ncol(cof_fwd[[length(mod_fwd)]])
Ftest[[length(mod_fwd)+1]]<-(rep(RSS_H0[[length(mod_fwd)+1]],length(RSSf[[length(mod_fwd)+1]]))/RSSf[[length(mod_fwd)+1]]-1)*df2[[length(mod_fwd)+1]]/df1
pval[[length(mod_fwd)+1]]<-pf(Ftest[[length(mod_fwd)+1]],df1,df2[[length(mod_fwd)+1]],lower.tail=FALSE)
rm(M,Y_t,cof_fwd_t,Res_H0,Q_,RSS)

##get max pval at each forward step
max_pval_fwd<-vector(mode="numeric",length=length(fwd_lm))
max_pval_fwd[1]<-0
for (i in 2:length(fwd_lm)) {max_pval_fwd[i]<-max(fwd_lm[[i]]$coef[(ncol(fix_cofs)+1):(ncol(fix_cofs)+ncol(cof_fwd[[i]])),4])}
rm(i)

##get the number of parameters & Loglikelihood from ML at each step
mod_fwd_LL<-list()
mod_fwd_LL[[1]]<-list(nfixed=ncol(cbind(fix_cofs,cof_fwd[[1]])),LL=emma.MLE(Y,cbind(fix_cofs,cof_fwd[[1]]),K_norm)$ML)
for (i in 2:length(cof_fwd)) {mod_fwd_LL[[i]]<-list(nfixed=ncol(cbind(fix_cofs,cof_fwd[[i]])),LL=emma.MLE(Y,cbind(fix_cofs,cof_fwd[[i]]),K_norm)$ML)}
rm(i)

cat('backward analysis','\n')

##BACKWARD (1st step == last fwd step)

dropcof_bwd<-list()
cof_bwd<-list()
mod_bwd <- list()
bwd_lm<-list()
herit_bwd<-list()

dropcof_bwd[[1]]<-'NA'
cof_bwd[[1]]<-as.matrix(cof_fwd[[length(mod_fwd)]][,!colnames(cof_fwd[[length(mod_fwd)]]) %in% dropcof_bwd[[1]]])
colnames(cof_bwd[[1]])<-colnames(cof_fwd[[length(mod_fwd)]])[!colnames(cof_fwd[[length(mod_fwd)]]) %in% dropcof_bwd[[1]]]
mod_bwd[[1]]<-emma.REMLE(Y,cbind(fix_cofs,cof_bwd[[1]]),K_norm)
herit_bwd[[1]]<-mod_bwd[[1]]$vg/(mod_bwd[[1]]$vg+mod_bwd[[1]]$ve)
M<-solve(chol(mod_bwd[[1]]$vg*K_norm+mod_bwd[[1]]$ve*diag(n)))
Y_t<-crossprod(M,Y)
cof_bwd_t<-crossprod(M,cbind(fix_cofs,cof_bwd[[1]]))
bwd_lm[[1]]<-summary(lm(Y_t~0+cof_bwd_t))

rm(M,Y_t,cof_bwd_t)


for (i in 2:length(mod_fwd)) {
dropcof_bwd[[i]]<-colnames(cof_bwd[[i-1]])[which(abs(bwd_lm[[i-1]]$coef[(ncol(fix_cofs)+1):nrow(bwd_lm[[i-1]]$coef),3])==min(abs(bwd_lm[[i-1]]$coef[(ncol(fix_cofs)+1):nrow(bwd_lm[[i-1]]$coef),3])))]
cof_bwd[[i]]<-as.matrix(cof_bwd[[i-1]][,!colnames(cof_bwd[[i-1]]) %in% dropcof_bwd[[i]]])
colnames(cof_bwd[[i]])<-colnames(cof_bwd[[i-1]])[!colnames(cof_bwd[[i-1]]) %in% dropcof_bwd[[i]]]
mod_bwd[[i]]<-emma.REMLE(Y,cbind(fix_cofs,cof_bwd[[i]]),K_norm)
herit_bwd[[i]]<-mod_bwd[[i]]$vg/(mod_bwd[[i]]$vg+mod_bwd[[i]]$ve)
M<-solve(chol(mod_bwd[[i]]$vg*K_norm+mod_bwd[[i]]$ve*diag(n)))
Y_t<-crossprod(M,Y)
cof_bwd_t<-crossprod(M,cbind(fix_cofs,cof_bwd[[i]]))
bwd_lm[[i]]<-summary(lm(Y_t~0+cof_bwd_t))
rm(M,Y_t,cof_bwd_t)}

rm(i)

##get max pval at each backward step
max_pval_bwd<-vector(mode="numeric",length=length(bwd_lm))
for (i in 1:(length(bwd_lm)-1)) {max_pval_bwd[i]<-max(bwd_lm[[i]]$coef[(ncol(fix_cofs)+1):(ncol(fix_cofs)+ncol(cof_bwd[[i]])),4])}
max_pval_bwd[length(bwd_lm)]<-0

##get the number of parameters & Loglikelihood from ML at each step
mod_bwd_LL<-list()
mod_bwd_LL[[1]]<-list(nfixed=ncol(cbind(fix_cofs,cof_bwd[[1]])),LL=emma.MLE(Y,cbind(fix_cofs,cof_bwd[[1]]),K_norm)$ML)
for (i in 2:length(cof_bwd)) {mod_bwd_LL[[i]]<-list(nfixed=ncol(cbind(fix_cofs,cof_bwd[[i]])),LL=emma.MLE(Y,cbind(fix_cofs,cof_bwd[[i]]),K_norm)$ML)}
rm(i)

cat('creating output','\n')

##Forward Table: Fwd + Bwd Tables
#Compute parameters for model criteria
BIC<-function(x){-2*x$LL+(x$nfixed+1)*log(n)}
extBIC<-function(x){BIC(x)+2*lchoose(m,x$nfixed-1)}

fwd_table<-data.frame(step=ncol(cof_fwd[[1]]),step_=paste('fwd',ncol(cof_fwd[[1]]),sep=''),cof=paste('+',addcof_fwd[[1]],sep=''),ncof=ncol(cof_fwd[[1]]),h2=herit_fwd[[1]]
	,maxpval=max_pval_fwd[1],BIC=BIC(mod_fwd_LL[[1]]),extBIC=extBIC(mod_fwd_LL[[1]]))
for (i in 2:(length(mod_fwd))) {fwd_table<-rbind(fwd_table,
	data.frame(step=ncol(cof_fwd[[i]]),step_=paste('fwd',ncol(cof_fwd[[i]]),sep=''),cof=paste('+',addcof_fwd[[i]],sep=''),ncof=ncol(cof_fwd[[i]]),h2=herit_fwd[[i]]
	,maxpval=max_pval_fwd[i],BIC=BIC(mod_fwd_LL[[i]]),extBIC=extBIC(mod_fwd_LL[[i]])))}

rm(i)

bwd_table<-data.frame(step=length(mod_fwd),step_=paste('bwd',0,sep=''),cof=paste('-',dropcof_bwd[[1]],sep=''),ncof=ncol(cof_bwd[[1]]),h2=herit_bwd[[1]]
	,maxpval=max_pval_bwd[1],BIC=BIC(mod_bwd_LL[[1]]),extBIC=extBIC(mod_bwd_LL[[1]]))
for (i in 2:(length(mod_bwd))) {bwd_table<-rbind(bwd_table,
	data.frame(step=length(mod_fwd)+i-1,step_=paste('bwd',i-1,sep=''),cof=paste('-',dropcof_bwd[[i]],sep=''),ncof=ncol(cof_bwd[[i]]),h2=herit_bwd[[i]]
	,maxpval=max_pval_bwd[i],BIC=BIC(mod_bwd_LL[[i]]),extBIC=extBIC(mod_bwd_LL[[i]])))}

rm(i,BIC,extBIC,max_pval_fwd,max_pval_bwd,dropcof_bwd)

fwdbwd_table<-rbind(fwd_table,bwd_table)

#RSS for plot

#null model only with intercept
null<-emma.REMLE(Y,as.matrix(Xo),K_norm)
M<-solve(chol(null$vg*K_norm+null$ve*diag(n)))
Y_t<-crossprod(M,Y)
Xo_t<-crossprod(M,as.matrix(Xo))
null_lm<-summary(lm(Y_t~0+Xo_t))
rm(null,M,Y_t,Xo_t)
RSS_null<-sum((Y-as.matrix(Xo)%*%null_lm$coef[,1])^2)

mod_fwd_RSS<-vector()
mod_fwd_RSS[1]<-sum((Y-cbind(fix_cofs,cof_fwd[[1]])%*%fwd_lm[[1]]$coef[,1])^2)
for (i in 2:length(mod_fwd)) {mod_fwd_RSS[i]<-sum((Y-cbind(fix_cofs,cof_fwd[[i]])%*%fwd_lm[[i]]$coef[,1])^2)}
mod_bwd_RSS<-vector()
mod_bwd_RSS[1]<-sum((Y-cbind(fix_cofs,cof_bwd[[1]])%*%bwd_lm[[1]]$coef[,1])^2)
for (i in 2:length(mod_bwd)) {mod_bwd_RSS[i]<-sum((Y-cbind(fix_cofs,cof_bwd[[i]])%*%bwd_lm[[i]]$coef[,1])^2)}

expl_RSS<-c(1-sapply(mod_fwd_RSS,function(x){x/RSS_null}),1-sapply(mod_bwd_RSS,function(x){x/RSS_null}))
fix_cofs_RSS<-rep(expl_RSS[1],length(expl_RSS))
cofs_RSS<-expl_RSS-fix_cofs_RSS
h2_RSS<-c(unlist(herit_fwd),unlist(herit_bwd))*(1-expl_RSS)
unexpl_RSS<-1-expl_RSS-h2_RSS
plot_RSS<-t(apply(cbind(fix_cofs_RSS,cofs_RSS,h2_RSS,unexpl_RSS),1,cumsum))

#GLS pvals at each step
pval_step<-list()
pval_step[[1]]<-list(out=data.frame('SNP'=names(pval[[2]]),'pval'=pval[[2]]),cof=addcof_fwd[[1]])
for (i in 2:(length(mod_fwd))) {pval_step[[i]]<-list('out'=rbind(data.frame('SNP'=colnames(cof_fwd[[i]]),'pval'=fwd_lm[[i]]$coef[(ncol(fix_cofs)+1):(ncol(fix_cofs)+ncol(cof_fwd[[i]])),4]),
	data.frame('SNP'=names(pval[[i+1]]),'pval'=pval[[i+1]])),cof=colnames(cof_fwd[[i]]))}

#GLS pvals for best models according to extBIC and mbonf

opt_extBIC<-fwdbwd_table[which(fwdbwd_table$extBIC==min(fwdbwd_table$extBIC))[1],]
opt_mbonf<-(fwdbwd_table[which(fwdbwd_table$maxpval<=0.05/m),])[which(fwdbwd_table[which(fwdbwd_table$maxpval<=0.05/m),]$ncof==max(fwdbwd_table[which(fwdbwd_table$maxpval<=0.05/m),]$ncof))[1],]
bestmodel_pvals<-function(model) {if(substr(model$step_,start=0,stop=3)=='fwd') {
		pval_step[[as.integer(substring(model$step_,first=4))+1]]} else if (substr(model$step_,start=0,stop=3)=='bwd') {
		cof<-cof_bwd[[as.integer(substring(model$step_,first=4))+1]]
		mixedmod<-emma.REMLE(Y,cbind(fix_cofs,cof),K_norm)
		M<-solve(chol(mixedmod$vg*K_norm+mixedmod$ve*diag(n)))
		Y_t<-crossprod(M,Y)
		cof_t<-crossprod(M,cbind(fix_cofs,cof))
		GLS_lm<-summary(lm(Y_t~0+cof_t))
		Res_H0<-GLS_lm$residuals
		Q_ <- qr.Q(qr(cof_t))
		RSS<-list()
		for (j in 1:(nbchunks-1)) {
		X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof)])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
		RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
		rm(X_t)}
		X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof)])[,((j)*round(m/nbchunks)+1):(m-ncol(cof))])
		RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
		rm(X_t,j)
		RSSf<-unlist(RSS)
		RSS_H0<-sum(Res_H0^2)
		df2<-n-df1-ncol(fix_cofs)-ncol(cof)
		Ftest<-(rep(RSS_H0,length(RSSf))/RSSf-1)*df2/df1
		pval<-pf(Ftest,df1,df2,lower.tail=FALSE)
		list(out=rbind(data.frame(SNP=colnames(cof),'pval'=GLS_lm$coef[(ncol(fix_cofs)+1):(ncol(fix_cofs)+ncol(cof)),4]),
		data.frame('SNP'=names(pval),'pval'=pval)),cof=colnames(cof))} else {cat('error \n')}}
opt_extBIC_out<-bestmodel_pvals(opt_extBIC)
opt_mbonf_out<-bestmodel_pvals(opt_mbonf)

list(step_table=fwdbwd_table,pval_step=pval_step,RSSout=plot_RSS,bonf_thresh=-log10(0.05/m),opt_extBIC=opt_extBIC_out,opt_mbonf=opt_mbonf_out)}

plot_step_table<-function(x,type){
	if (type=='h2') {plot(x$step_table$step,x$step_table$h2,type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='h2')
	abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
    else if (type=='maxpval'){plot(x$step_table$step,-log10(x$step_table$maxpval),type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='-log10(max_Pval)')
	abline(h=x$bonf_thresh,lty=2)
	abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
    else if (type=='BIC'){plot(x$step_table$step,x$step_table$BIC,type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='BIC')
	abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
    else if (type=='extBIC'){plot(x$step_table$step,x$step_table$extBIC,type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='EBIC')
	abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
	else {cat('error! \n argument type must be one of h2, maxpval, BIC, extBIC')}}

plot_step_RSS<-function(x){
	op<-par(mar=c(5, 5, 2, 2))
	plot(0,0,xlim=c(0,nrow(x$RSSout)-1),ylim=c(0,1),xlab='step',ylab='%var',col=0)
	polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,4],0,0), col='brown1', border=0)
	polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,3],0,0), col='forestgreen', border=0)
	polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,2],0,0), col='dodgerblue4', border=0)
	polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,1],0,0), col='grey', border=0)
	abline(v=(nrow(x$step_table)/2-0.5),lty=2)
	par(op)}

plot_GWAS<-function(x) {
	output_<-x$out[order(x$out$Pos),]
	output_ok<-output_[order(output_$Chr),]
	
	maxpos<-c(0,cumsum(as.numeric(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x+max(cumsum(as.numeric(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x)))/200)))
	plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$Chr))/2))
#	plot_col<-c('blue','darkgreen','red','cyan','purple')
	size<-aggregate(output_ok$Pos,list(output_ok$Chr),length)$x
	a<-rep(maxpos[1],size[1])
	b<-rep(plot_col[1],size[1])
		for (i in 2:max(unique(output_ok$Chr))){
	a<-c(a,rep(maxpos[i],size[i]))
	b<-c(b,rep(plot_col[i],size[i]))}

	output_ok$xpos<-output_ok$Pos+a
	output_ok$col<-b
	output_ok$col[output_ok$SNP %in% x$cof]<-'red'
	
	d<-(aggregate(output_ok$xpos,list(output_ok$Chr),min)$x+aggregate(output_ok$xpos,list(output_ok$Chr),max)$x)/2

	plot(output_ok$xpos,-log10(output_ok$pval),col=output_ok$col,pch=20,ylab='-log10(pval)',xaxt='n',xlab='chromosome')
	axis(1,tick=FALSE,at=d,labels=c(1:max(unique(output_ok$Chr))))
	abline(h=x$bonf_thresh,lty=3,col='black')}

plot_region<-function(x,chrom,pos1,pos2){
	region<-subset(x$out,Chr==chrom & Pos>=pos1 & Pos <=pos2)
	region$col<- if (chrom %% 2 == 0) {'gray60'} else {'gray10'}
	region$col[which(region$SNP %in% x$cof)]<-'red'
	plot(region$Pos,-log10(region$pval),type='p',pch=20,main=paste('chromosome',chrom,sep=''),xlab='position (bp)',ylab='-log10(pval)',col=region$col,xlim=c(pos1,pos2))
	abline(h=x$bonf_thresh,lty=3,col='black')}


plot_fwd_GWAS<-function(x,step,snp_info,pval_filt) {
	stopifnot(step<=length(x$pval_step))
	output<-list(out=subset(merge(snp_info,x$pval_step[[step]]$out,by='SNP'),pval<=pval_filt),cof=x$pval_step[[step]]$cof,bonf_thresh=x$bonf_thresh)
	plot_GWAS(output)}

plot_fwd_region<-function(x,step,snp_info,pval_filt,chrom,pos1,pos2) {
	stopifnot(step<=length(x$pval_step))
	output<-list(out=subset(merge(snp_info,x$pval_step[[step]]$out,by='SNP'),pval<=pval_filt),cof=x$pval_step[[step]]$cof,bonf_thresh=x$bonf_thresh)
	plot_region(output,chrom,pos1,pos2)}


plot_opt_GWAS<-function(x,opt,snp_info,pval_filt) {
	if (opt=='extBIC') {output<-list(out=subset(merge(snp_info,x$opt_extBIC$out,by='SNP'),pval<=pval_filt),cof=x$opt_extBIC$cof,bonf_thresh=x$bonf_thresh)
		plot_GWAS(output)}
	else if (opt=='mbonf') {output<-list(out=subset(merge(snp_info,x$opt_mbonf$out,by='SNP'),pval<=pval_filt),cof=x$opt_mbonf$cof,bonf_thresh=x$bonf_thresh)
		plot_GWAS(output)}
	else {cat('error! \n opt must be extBIC or mbonf')}}

plot_opt_region<-function(x,opt,snp_info,pval_filt,chrom,pos1,pos2) {
	if (opt=='extBIC') {output<-list(out=subset(merge(snp_info,x$opt_extBIC$out,by='SNP'),pval<=pval_filt),cof=x$opt_extBIC$cof,bonf_thresh=x$bonf_thresh)
		plot_region(output,chrom,pos1,pos2)}
	else if (opt=='mbonf') {output<-list(out=subset(merge(snp_info,x$opt_mbonf$out,by='SNP'),pval<=pval_filt),cof=x$opt_mbonf$cof,bonf_thresh=x$bonf_thresh)
		plot_region(output,chrom,pos1,pos2)}
	else {cat('error! \n opt must be extBIC or mbonf')}}


qqplot_fwd_GWAS<-function(x,nsteps){
	stopifnot(nsteps<=length(x$pval_step))
	e<--log10(ppoints(nrow(x$pval_step[[1]]$out)))
	ostep<-list()
	ostep[[1]]<--log10(sort(x$pval_step[[1]]$out$pval))
	for (i in 2:nsteps) {ostep[[i]]<--log10(sort(x$pval_step[[i]]$out$pval))}

	maxp<-ceiling(max(unlist(ostep)))

	plot(e,ostep[[1]],type='b',pch=20,cex=0.8,col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp))
	abline(0,1,col="dark grey")

	for (i in 2:nsteps) {
	par(new=T)
	plot(e,ostep[[i]],type='b',pch=20,cex=0.8,col=i,axes='F',xlab='',ylab='',xlim=c(0,max(e)+1),ylim=c(0,maxp))}
	legend(0,maxp,lty=1,pch=20,col=c(1:length(ostep)),paste(c(0:(length(ostep)-1)),'cof',sep=' '))
}

qqplot_opt_GWAS<-function(x,opt){
	if (opt=='extBIC') {
		e<--log10(ppoints(nrow(x$opt_extBIC$out)))
		o<--log10(sort(x$opt_extBIC$out$pval))
		maxp<-ceiling(max(o))
		plot(e,o,type='b',pch=20,cex=0.8,col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),main=paste('optimal model according to extBIC'))
		abline(0,1,col="dark grey")}
	else if (opt=='mbonf') {
		e<--log10(ppoints(nrow(x$opt_mbonf$out)))
		o<--log10(sort(x$opt_mbonf$out$pval))
		maxp<-ceiling(max(o))
		plot(e,o,type='b',pch=20,cex=0.8,col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),main=paste('optimal model according to mbonf'))
		abline(0,1,col="dark grey")}
	else {cat('error! \n opt must be extBIC or mbonf')}}


