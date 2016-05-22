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

##' MLMM
##'
##' MLMM
##' @param Y phenotypes, a vector of length m, with names(Y)=individual names
##' @param X genotypes, a n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names
##' @param K kinship, a n by n matrix, with rownames(K)=colnames(K)=individual names
##' @param nbchunks an integer defining the number of chunks of X to run the analysis, allows to decrease the memory usage ==> minimum=2, increase it if you do not have enough memory
##' @param maxsteps maximum number of steps desired in the forward approach. The forward approach breaks automatically once the pseudo-heritability is close to 0, however to avoid doing too many steps in case the pseudo-heritability does not reach a value close to 0, this parameter is also used. It's value must be specified as an integer >= 3
##' @param thresh threshold
##' @return results
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
mlmm <- function(Y,X,K,nbchunks,maxsteps,thresh = NULL) {

  n<-length(Y)
  m<-ncol(X)

  stopifnot(ncol(K) == n)
  stopifnot(nrow(K) == n)
  stopifnot(nrow(X) == n)
  stopifnot(nbchunks >= 2)
  stopifnot(maxsteps >= 3)

  ##INTERCEPT

  Xo<-rep(1,n)

  ##K MATRIX NORMALISATION

  K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
  rm(K)

  ##step 0 : NULL MODEL
  cof_fwd<-list()
  cof_fwd[[1]]<-as.matrix(Xo)
  colnames(cof_fwd[[1]])<-'Xo'

  mod_fwd<-list()
  mod_fwd[[1]]<-emma::emma.REMLE(Y,cof_fwd[[1]],K_norm)

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

  ##step 1 : EMMAX

  M<-solve(chol(mod_fwd[[1]]$vg*K_norm+mod_fwd[[1]]$ve*diag(n)))
  Y_t<-crossprod(M,Y)
  cof_fwd_t<-crossprod(M,cof_fwd[[1]])
  fwd_lm[[1]]<-summary(lm(Y_t~0+cof_fwd_t))
  Res_H0<-fwd_lm[[1]]$residuals
  Q_<-qr.Q(qr(cof_fwd_t))

  RSS<-list()
  for (j in 1:(nbchunks-1)) {
    X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[1]])])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
    rm(X_t)}
  X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[1]])])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof_fwd[[1]])-1))])
  RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
  rm(X_t,j)

  RSSf[[2]]<-unlist(RSS)
  RSS_H0[[2]]<-sum(Res_H0^2)
  df2[[2]]<-n-df1-ncol(cof_fwd[[1]])
  Ftest[[2]]<-(rep(RSS_H0[[2]],length(RSSf[[2]]))/RSSf[[2]]-1)*df2[[2]]/df1
  pval[[2]]<-pf(Ftest[[2]],df1,df2[[2]],lower.tail=FALSE)

  cof_fwd[[2]]<-cbind(cof_fwd[[1]],X[,colnames(X) %in% names(which(RSSf[[2]]==min(RSSf[[2]]))[1])])
  colnames(cof_fwd[[2]])<-c(colnames(cof_fwd[[1]]),names(which(RSSf[[2]]==min(RSSf[[2]]))[1]))
  mod_fwd[[2]]<-emma::emma.REMLE(Y,cof_fwd[[2]],K_norm)
  herit_fwd[[2]]<-mod_fwd[[2]]$vg/(mod_fwd[[2]]$vg+mod_fwd[[2]]$ve)
  rm(M,Y_t,cof_fwd_t,Res_H0,Q_,RSS)

  cat('step 1 done! pseudo-h=',round(herit_fwd[[2]],3),'\n')

  ##FORWARD

  for (i in 3:(maxsteps)) {
    if (herit_fwd[[i-2]] < 0.01){
      break
    } else {

      M<-solve(chol(mod_fwd[[i-1]]$vg*K_norm+mod_fwd[[i-1]]$ve*diag(n)))
      Y_t<-crossprod(M,Y)
      cof_fwd_t<-crossprod(M,cof_fwd[[i-1]])
      fwd_lm[[i-1]]<-summary(lm(Y_t~0+cof_fwd_t))
      Res_H0<-fwd_lm[[i-1]]$residuals
      Q_ <- qr.Q(qr(cof_fwd_t))

      RSS<-list()
      for (j in 1:(nbchunks-1)) {
        X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[i-1]])])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
        RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
        rm(X_t)}
      X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[i-1]])])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof_fwd[[i-1]])-1))])
      RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
      rm(X_t,j)

      RSSf[[i]]<-unlist(RSS)
      RSS_H0[[i]]<-sum(Res_H0^2)
      df2[[i]]<-n-df1-ncol(cof_fwd[[i-1]])
      Ftest[[i]]<-(rep(RSS_H0[[i]],length(RSSf[[i]]))/RSSf[[i]]-1)*df2[[i]]/df1
      pval[[i]]<-pf(Ftest[[i]],df1,df2[[i]],lower.tail=FALSE)

      cof_fwd[[i]]<-cbind(cof_fwd[[i-1]],X[,colnames(X) %in% names(which(RSSf[[i]]==min(RSSf[[i]]))[1])])
      colnames(cof_fwd[[i]])<-c(colnames(cof_fwd[[i-1]]),names(which(RSSf[[i]]==min(RSSf[[i]]))[1]))
      mod_fwd[[i]]<-emma::emma.REMLE(Y,cof_fwd[[i]],K_norm)
      herit_fwd[[i]]<-mod_fwd[[i]]$vg/(mod_fwd[[i]]$vg+mod_fwd[[i]]$ve)
      rm(M,Y_t,cof_fwd_t,Res_H0,Q_,RSS)}
    cat('step ',i-1,' done! pseudo-h=',round(herit_fwd[[i]],3),'\n')}
  rm(i)

  ##gls at last forward step
  M<-solve(chol(mod_fwd[[length(mod_fwd)]]$vg*K_norm+mod_fwd[[length(mod_fwd)]]$ve*diag(n)))
  Y_t<-crossprod(M,Y)
  cof_fwd_t<-crossprod(M,cof_fwd[[length(mod_fwd)]])
  fwd_lm[[length(mod_fwd)]]<-summary(lm(Y_t~0+cof_fwd_t))

  Res_H0<-fwd_lm[[length(mod_fwd)]]$residuals
  Q_ <- qr.Q(qr(cof_fwd_t))

  RSS<-list()
  for (j in 1:(nbchunks-1)) {
    X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[length(mod_fwd)]])])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
    rm(X_t)}
  X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof_fwd[[length(mod_fwd)]])])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof_fwd[[length(mod_fwd)]])-1))])
  RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
  rm(X_t,j)

  RSSf[[length(mod_fwd)+1]]<-unlist(RSS)
  RSS_H0[[length(mod_fwd)+1]]<-sum(Res_H0^2)
  df2[[length(mod_fwd)+1]]<-n-df1-ncol(cof_fwd[[length(mod_fwd)]])
  Ftest[[length(mod_fwd)+1]]<-(rep(RSS_H0[[length(mod_fwd)+1]],length(RSSf[[length(mod_fwd)+1]]))/RSSf[[length(mod_fwd)+1]]-1)*df2[[length(mod_fwd)+1]]/df1
  pval[[length(mod_fwd)+1]]<-pf(Ftest[[length(mod_fwd)+1]],df1,df2[[length(mod_fwd)+1]],lower.tail=FALSE)
  rm(M,Y_t,cof_fwd_t,Res_H0,Q_,RSS)

  ##get max pval at each forward step
  max_pval_fwd<-vector(mode="numeric",length=length(fwd_lm))
  max_pval_fwd[1]<-0
  for (i in 2:length(fwd_lm)) {max_pval_fwd[i]<-max(fwd_lm[[i]]$coef[2:i,4])}
  rm(i)

  ##get the number of parameters & Loglikelihood from ML at each step
  mod_fwd_LL<-list()
  mod_fwd_LL[[1]]<-list(nfixed=ncol(cof_fwd[[1]]),LL=emma::emma.MLE(Y,cof_fwd[[1]],K_norm)$ML)
  for (i in 2:length(cof_fwd)) {mod_fwd_LL[[i]]<-list(nfixed=ncol(cof_fwd[[i]]),LL=emma::emma.MLE(Y,cof_fwd[[i]],K_norm)$ML)}
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
  mod_bwd[[1]]<-emma::emma.REMLE(Y,cof_bwd[[1]],K_norm)
  herit_bwd[[1]]<-mod_bwd[[1]]$vg/(mod_bwd[[1]]$vg+mod_bwd[[1]]$ve)
  M<-solve(chol(mod_bwd[[1]]$vg*K_norm+mod_bwd[[1]]$ve*diag(n)))
  Y_t<-crossprod(M,Y)
  cof_bwd_t<-crossprod(M,cof_bwd[[1]])
  bwd_lm[[1]]<-summary(lm(Y_t~0+cof_bwd_t))

  rm(M,Y_t,cof_bwd_t)

  for (i in 2:length(mod_fwd)) {
    dropcof_bwd[[i]]<-(colnames(cof_bwd[[i-1]])[2:ncol(cof_bwd[[i-1]])])[which(abs(bwd_lm[[i-1]]$coef[2:nrow(bwd_lm[[i-1]]$coef),3])==min(abs(bwd_lm[[i-1]]$coef[2:nrow(bwd_lm[[i-1]]$coef),3])))]
    cof_bwd[[i]]<-as.matrix(cof_bwd[[i-1]][,!colnames(cof_bwd[[i-1]]) %in% dropcof_bwd[[i]]])
    colnames(cof_bwd[[i]])<-colnames(cof_bwd[[i-1]])[!colnames(cof_bwd[[i-1]]) %in% dropcof_bwd[[i]]]
    mod_bwd[[i]]<-emma::emma.REMLE(Y,cof_bwd[[i]],K_norm)
    herit_bwd[[i]]<-mod_bwd[[i]]$vg/(mod_bwd[[i]]$vg+mod_bwd[[i]]$ve)
    M<-solve(chol(mod_bwd[[i]]$vg*K_norm+mod_bwd[[i]]$ve*diag(n)))
    Y_t<-crossprod(M,Y)
    cof_bwd_t<-crossprod(M,cof_bwd[[i]])
    bwd_lm[[i]]<-summary(lm(Y_t~0+cof_bwd_t))
    rm(M,Y_t,cof_bwd_t)}

  rm(i)

  ##get max pval at each backward step
  max_pval_bwd<-vector(mode="numeric",length=length(bwd_lm))
  for (i in 1:(length(bwd_lm)-1)) {max_pval_bwd[i]<-max(bwd_lm[[i]]$coef[2:(length(bwd_lm)+1-i),4])}
  max_pval_bwd[length(bwd_lm)]<-0

  ##get the number of parameters & Loglikelihood from ML at each step
  mod_bwd_LL<-list()
  mod_bwd_LL[[1]]<-list(nfixed=ncol(cof_bwd[[1]]),LL=emma::emma.MLE(Y,cof_bwd[[1]],K_norm)$ML)
  for (i in 2:length(cof_bwd)) {mod_bwd_LL[[i]]<-list(nfixed=ncol(cof_bwd[[i]]),LL=emma::emma.MLE(Y,cof_bwd[[i]],K_norm)$ML)}
  rm(i)

  cat('creating output','\n')

  ##Forward Table: Fwd + Bwd Tables
  ##Compute parameters for model criteria
  BIC<-function(x){-2*x$LL+(x$nfixed+1)*log(n)}
  extBIC<-function(x){BIC(x)+2*lchoose(m,x$nfixed-1)}

  fwd_table<-data.frame(step=ncol(cof_fwd[[1]])-1,step_=paste('fwd',ncol(cof_fwd[[1]])-1,sep=''),cof='NA',ncof=ncol(cof_fwd[[1]])-1,h2=herit_fwd[[1]]
                       ,maxpval=max_pval_fwd[1],BIC=BIC(mod_fwd_LL[[1]]),extBIC=extBIC(mod_fwd_LL[[1]]))
  for (i in 2:(length(mod_fwd))) {fwd_table<-rbind(fwd_table,
                                                   data.frame(step=ncol(cof_fwd[[i]])-1,step_=paste('fwd',ncol(cof_fwd[[i]])-1,sep=''),cof=paste('+',colnames(cof_fwd[[i]])[i],sep=''),ncof=ncol(cof_fwd[[i]])-1,h2=herit_fwd[[i]]
                                                             ,maxpval=max_pval_fwd[i],BIC=BIC(mod_fwd_LL[[i]]),extBIC=extBIC(mod_fwd_LL[[i]])))}

  rm(i)

  bwd_table<-data.frame(step=length(mod_fwd),step_=paste('bwd',0,sep=''),cof=paste('-',dropcof_bwd[[1]],sep=''),ncof=ncol(cof_bwd[[1]])-1,h2=herit_bwd[[1]]
                       ,maxpval=max_pval_bwd[1],BIC=BIC(mod_bwd_LL[[1]]),extBIC=extBIC(mod_bwd_LL[[1]]))
  for (i in 2:(length(mod_bwd))) {bwd_table<-rbind(bwd_table,
                                                   data.frame(step=length(mod_fwd)+i-1,step_=paste('bwd',i-1,sep=''),cof=paste('-',dropcof_bwd[[i]],sep=''),ncof=ncol(cof_bwd[[i]])-1,h2=herit_bwd[[i]]
                                                             ,maxpval=max_pval_bwd[i],BIC=BIC(mod_bwd_LL[[i]]),extBIC=extBIC(mod_bwd_LL[[i]])))}

  rm(i,BIC,extBIC,max_pval_fwd,max_pval_bwd,dropcof_bwd)

  fwdbwd_table<-rbind(fwd_table,bwd_table)

  ##RSS for plot
  mod_fwd_RSS<-vector()
  mod_fwd_RSS[1]<-sum((Y-cof_fwd[[1]]%*%fwd_lm[[1]]$coef[,1])^2)
  for (i in 2:length(mod_fwd)) {mod_fwd_RSS[i]<-sum((Y-cof_fwd[[i]]%*%fwd_lm[[i]]$coef[,1])^2)}
  mod_bwd_RSS<-vector()
  mod_bwd_RSS[1]<-sum((Y-cof_bwd[[1]]%*%bwd_lm[[1]]$coef[,1])^2)
  for (i in 2:length(mod_bwd)) {mod_bwd_RSS[i]<-sum((Y-cof_bwd[[i]]%*%bwd_lm[[i]]$coef[,1])^2)}
  expl_RSS<-c(1-sapply(mod_fwd_RSS,function(x){x/mod_fwd_RSS[1]}),1-sapply(mod_bwd_RSS,function(x){x/mod_bwd_RSS[length(mod_bwd_RSS)]}))
  h2_RSS<-c(unlist(herit_fwd),unlist(herit_bwd))*(1-expl_RSS)
  unexpl_RSS<-1-expl_RSS-h2_RSS
  plot_RSS<-t(apply(cbind(expl_RSS,h2_RSS,unexpl_RSS),1,cumsum))

  ##GLS pvals at each step
  pval_step<-list()
  pval_step[[1]]<-list(out=data.frame("SNP"=colnames(X),"pval"=pval[[2]]),"cof"=NA, "coef"=fwd_lm[[1]]$coef)
  for (i in 2:(length(mod_fwd))) {pval_step[[i]]<-list(out=rbind(data.frame(SNP=colnames(cof_fwd[[i]])[-1],'pval'=fwd_lm[[i]]$coef[2:i,4]),
                                                                 data.frame(SNP=colnames(X)[-which(colnames(X) %in% colnames(cof_fwd[[i]]))],'pval'=pval[[i+1]])),"cof"=colnames(cof_fwd[[i]])[-1], "coef"=fwd_lm[[i]]$coef)}

  ##GLS pvals for best models according to extBIC and mbonf

  opt_extBIC<-fwdbwd_table[which(fwdbwd_table$extBIC==min(fwdbwd_table$extBIC))[1],]
  opt_mbonf<-(fwdbwd_table[which(fwdbwd_table$maxpval<=0.05/m),])[which(fwdbwd_table[which(fwdbwd_table$maxpval<=0.05/m),]$ncof==max(fwdbwd_table[which(fwdbwd_table$maxpval<=0.05/m),]$ncof))[1],]
  if(! is.null(thresh)){
    opt_thresh<-(fwdbwd_table[which(fwdbwd_table$maxpval<=thresh),])[which(fwdbwd_table[which(fwdbwd_table$maxpval<=thresh),]$ncof==max(fwdbwd_table[which(fwdbwd_table$maxpval<=thresh),]$ncof))[1],]
  }
  bestmodel_pvals<-function(model) {
    if(substr(model$step_,start=0,stop=3)=='fwd') {
      pval_step[[as.integer(substring(model$step_,first=4))+1]]
    } else if (substr(model$step_,start=0,stop=3)=='bwd') {
      cof<-cof_bwd[[as.integer(substring(model$step_,first=4))+1]]
      mixedmod<-emma::emma.REMLE(Y,cof,K_norm)
      M<-solve(chol(mixedmod$vg*K_norm+mixedmod$ve*diag(n)))
      Y_t<-crossprod(M,Y)
      cof_t<-crossprod(M,cof)
      GLS_lm<-summary(lm(Y_t~0+cof_t))
      Res_H0<-GLS_lm$residuals
      Q_ <- qr.Q(qr(cof_t))
      RSS<-list()
      for (j in 1:(nbchunks-1)) {
        X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof)])[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
        RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
        rm(X_t)}
      X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),(X[,!colnames(X) %in% colnames(cof)])[,((j)*round(m/nbchunks)+1):(m-(ncol(cof)-1))])
      RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
      rm(X_t,j)
      RSSf<-unlist(RSS)
      RSS_H0<-sum(Res_H0^2)
      df2<-n-df1-ncol(cof)
      Ftest<-(rep(RSS_H0,length(RSSf))/RSSf-1)*df2/df1
      pval<-pf(Ftest,df1,df2,lower.tail=FALSE)
      list('out'=rbind(data.frame(SNP=colnames(cof)[-1],'pval'=GLS_lm$coef[2:(ncol(cof)),4]),
                       data.frame('SNP'=colnames(X)[-which(colnames(X) %in% colnames(cof))],'pval'=pval)),
           'cof'=colnames(cof)[-1],
           'coef'=GLS_lm$coef)} else {cat('error \n')}}
  opt_extBIC_out<-bestmodel_pvals(opt_extBIC)
  opt_mbonf_out<-bestmodel_pvals(opt_mbonf)
  if(! is.null(thresh)){
    opt_thresh_out<-bestmodel_pvals(opt_thresh)
  }
  output <- list(step_table=fwdbwd_table,pval_step=pval_step,RSSout=plot_RSS,bonf_thresh=-log10(0.05/m),opt_extBIC=opt_extBIC_out,opt_mbonf=opt_mbonf_out)
  if(! is.null(thresh)){
    output$thresh <- -log10(thresh)
    output$opt_thresh <- opt_thresh_out
  }
  return(output)
}
