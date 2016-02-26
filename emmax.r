##############################################################################################################################################
###EMMAX
###SET OF FUNCTIONS TO RUN GWAS CORRECTING FOR POPULATION STRUCTURE WITH EMMAX (Kang et al. 2010, NatGen 42:348-354)
#######
#
##note: require EMMA
#library(emma)
#source('emma.r')
#
##REQUIRED DATA & FORMAT
#
#PHENOTYPE - Y: a vector of length n with names(Y)=ecotype names
#GENOTYPE - X: a n by m matrix, where n=number of ecotypes, m=number of markers, with rownames(X)=ecotype names, and colnames(X)=SNP names
#KINSHIP - K: a n by n matrix, with rownames(K)=colnames(K)=ecotype names
#each of these data being sorted in the same way, according to the ecotype name
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
#source('path/fwd_emmax.r')
#
###EMMAX SCAN
#mygwas<-emmax(Y,X,K,nbchunks)
#X,Y,K as described above
#nbchunks: an integer defining the number of chunks of X to run the analysis, allows to decrease the memory usage ==> minimum=2, increase it if you do not have enough memory 
#
###RESULTS
#
##FUNCTION OUTPUT
#A LIST:
#	$output: a data.frame with the F statitics, pvals and R2 for each SNP tested
#	$bonf_thresh: pval threshold according to a bonferonni correction for an alpha of 0.05
#
##PLOTS
#
#GWAS MANHATTAN PLOT
#plot_GWAS(mygwas,snp_info,pval_filt)
#snp_info as described above
#pval_filt=a p-value threshold for filtering the output, only p-vals below this threshold will be displayed in the plot
#
#GWAS MANHATTAN PLOT ZOOMED IN A REGION OF INTEREST
#plot_region(mygwas,snp_info,chrom,pos1,pos2)
#step, snp_info as described above
#chrom=on which chromosome is the region of interest
#pos1, pos2=delimitations of the region of interest in the same unit as Pos in snp_info
#
#p-values QQplot
#qqplot_GWAS(mygwas)
##############################################################################################################################################

emmax<-function(Y,X,K,nbchunks) {

n<-length(Y)
m<-ncol(X)

stopifnot(ncol(K) == n)
stopifnot(nrow(K) == n)
stopifnot(nrow(X) == n)
stopifnot(nbchunks >= 2)

#INTERCEPT

Xo<-rep(1,n)

#K MATRIX NORMALISATION

K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
rm(K)

#NULL MODEL

null<-emma.REMLE(Y,as.matrix(Xo),K_norm)

pseudoh<-null$vg/(null$vg+null$ve)

cat('null model done! pseudo-h =',round(pseudoh,3),'\n')

#EMMAX

M<-solve(chol(null$vg*K_norm+null$ve*diag(n)))
Y_t<-crossprod(M,Y)
Xo_t<-crossprod(M,Xo)

RSS<-list()
for (j in 1:(nbchunks-1)) {
X_t<-crossprod(M,X[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(cbind(Xo_t,x),Y_t,intercept = FALSE)$residuals^2)})
rm(X_t)}
X_t<-crossprod(M,X[,((j)*round(m/nbchunks)+1):(m)])
RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(cbind(Xo_t,x),Y_t,intercept = FALSE)$residuals^2)})
rm(X_t,j)

RSSf<-unlist(RSS)
RSS_H0<-rep(sum(lsfit(Xo_t,Y_t,intercept = FALSE)$residuals^2),m)
df1<-1
df2<-n-df1-1
R2<-1-1/(RSS_H0/RSSf)
F<-(RSS_H0/RSSf-1)*df2/df1
pval<-pf(F,df1,df2,lower.tail=FALSE)

cat('EMMAX scan done! \n')

cat('creating output','\n')

list(output=data.frame(SNP=colnames(X),'F'=F,'pval'=pval,'Rsq'=R2),bonf_thresh=-log10(0.05/m))}

linreg<-function(Y,X,nbchunks) {

n<-length(Y)
m<-ncol(X)

stopifnot(nrow(X) == n)
stopifnot(nbchunks >= 2)

#INTERCEPT

Xo<-rep(1,n)

RSS<-list()
for (j in 1:(nbchunks-1)) {RSS[[j]]<-apply(X[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))],2,function(x){sum(lsfit(cbind(Xo,x),Y,intercept=FALSE)$residuals^2)})}
RSS[[nbchunks]]<-apply(X[,((j)*round(m/nbchunks)+1):(m)],2,function(x){sum(lsfit(cbind(Xo,x),Y,intercept=FALSE)$residuals^2)})
rm(j)

RSSf<-unlist(RSS)
RSS_H0<-rep(sum(lsfit(Xo,Y,intercept=FALSE)$residuals^2),m)
df1<-1
df2<-n-df1-1
R2<-1-1/(RSS_H0/RSSf)
F<-(RSS_H0/RSSf-1)*df2/df1
pval<-pf(F,df1,df2,lower.tail=FALSE)

cat('linreg scan done! \n')

cat('creating output','\n')

list(output=data.frame(SNP=colnames(X),'F'=F,'pval'=pval,'Rsq'=R2),bonf_thresh=-log10(0.05/m))}


plot_GWAS<-function(x,snp_info,pval_filt) {

output<-subset(merge(snp_info,x$output,by='SNP'),pval<=pval_filt)
output_<-output[order(output$Pos),]
output_ok<-output_[order(output_$Chr),]

maxpos<-c(0,cumsum(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x+max(cumsum(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x))/100))
plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$Chr))/2))
#plot_col<-c('blue','darkgreen','red','cyan','purple')
size<-aggregate(output_ok$Pos,list(output_ok$Chr),length)$x

a<-rep(maxpos[1],size[1])
b<-rep(plot_col[1],size[1])
if (length(unique(output_ok$Chr))>1){
	for (i in 2:length(unique(output_ok$Chr))){
a<-c(a,rep(maxpos[i],size[i]))
b<-c(b,rep(plot_col[i],size[i]))}}

output_ok$xpos<-output_ok$Pos+a
output_ok$col<-b

d<-(aggregate(output_ok$xpos,list(output_ok$Chr),min)$x+aggregate(output_ok$xpos,list(output_ok$Chr),max)$x)/2

plot(output_ok$xpos,-log10(output_ok$pval),col=output_ok$col,pch=20,ylab='-log10(pval)',xaxt='n',xlab='chromosome')
axis(1,tick=FALSE,at=d,labels=unique(output_ok$Chr))
abline(h=x$bonf_thresh,lty=3,col='black')}


plot_region<-function(x,snp_info,chrom,pos1,pos2){

output<-merge(snp_info,x$output,by='SNP')
region<-subset(output,Chr==chrom & Pos>=pos1 & Pos <=pos2)

plot(region$Pos,-log10(region$pval),type='p',pch=20,main=paste('chromosome',chrom,sep=''),xlab='position (bp)',ylab='-log10(pval)',col='gray40',xlim=c(pos1,pos2))
abline(h=x$bonf_thresh,lty=3,col='black')}

qqplot_GWAS<-function(x){
e<--log10(ppoints(nrow(x$output)))
o<--log10(sort(x$output$pval))

plot(e,o,type='b',pch=20,cex=0.8,col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))))
abline(0,1,col="dark grey")}
