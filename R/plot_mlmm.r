##' Plot
##'
##' Plot
##' @param x x
##' @param type type
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, such as \code{main}
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_step_table<-function(x,type,...){
  if (type=='h2') {graphics::plot(x$step_table$step,x$step_table$h2,type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='h2')
    graphics::abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
  else if (type=='maxpval'){graphics::plot(x$step_table$step,-log10(x$step_table$maxpval),type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='-log10(max_Pval)')
    graphics::abline(h=x$bonf_thresh,lty=2)
    if(! is.null(x$thresh)){graphics::abline(h=x$thresh,lty=2,col=2)}
    graphics::abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
  else if (type=='BIC'){graphics::plot(x$step_table$step,x$step_table$BIC,type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='BIC')
    graphics::abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
  else if (type=='extBIC'){graphics::plot(x$step_table$step,x$step_table$extBIC,type='b',lty=2,pch=20,col='darkblue',xlab='step',ylab='EBIC')
    graphics::abline(v=(nrow(x$step_table)/2-0.5),lty=2)}
  else {cat('error! \n argument type must be one of h2, maxpval, BIC, extBIC')}}

##' Plot
##'
##' Plot
##' @param x x
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, such as \code{main}
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_step_RSS<-function(x,...){
  op<-graphics::par(mar=c(5, 5, 2, 2))
  graphics::plot(0,0,xlim=c(0,nrow(x$RSSout)-1),ylim=c(0,1),xlab='step',ylab='%var',col=0,...)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,3],0,0), col='brown1', border=0)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,2],0,0), col='forestgreen', border=0)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,1],0,0), col='dodgerblue4', border=0)
  graphics::abline(v=(nrow(x$step_table)/2-0.5),lty=2)
  graphics::par(op)}

##' Plot
##'
##' Plot
##' @param x x
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, such as \code{main}
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_step_RSS_cof<-function(x,...){
  op<-graphics::par(mar=c(5, 5, 2, 2))
  graphics::plot(0,0,xlim=c(0,nrow(x$RSSout)-1),ylim=c(0,1),xlab='step',ylab='%var',col=0,...)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,4],0,0), col='brown1', border=0)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,3],0,0), col='forestgreen', border=0)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,2],0,0), col='dodgerblue4', border=0)
  graphics::polygon(c(0:(nrow(x$RSSout)-1),(nrow(x$RSSout)-1),0), c(x$RSSout[,1],0,0), col='grey', border=0)
  graphics::abline(v=(nrow(x$step_table)/2-0.5),lty=2)
  graphics::par(op)}

##' Plot
##'
##' Plot
##' @param x x
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, such as \code{main}
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_GWAS<-function(x,...) {
  output_<-x$out[order(x$out$Pos),]
  output_ok<-output_[order(output_$Chr),]
  maxpos<-c(0,cumsum(as.numeric(stats::aggregate(output_ok$Pos,list(output_ok$Chr),max)$x+max(cumsum(as.numeric(stats::aggregate(output_ok$Pos,list(output_ok$Chr),max)$x)))/200)))
  plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$Chr))/2))
                                        #	plot_col<-c('blue','darkgreen','red','cyan','purple')
  size<-stats::aggregate(output_ok$Pos,list(output_ok$Chr),length)$x
  a<-rep(maxpos[1],size[1])
  b<-rep(plot_col[1],size[1])
  if (length(unique(output_ok$Chr))>1){
    for (i in 2:length(unique(output_ok$Chr))){
      a<-c(a,rep(maxpos[i],size[i]))
      b<-c(b,rep(plot_col[i],size[i]))}}
  output_ok$xpos<-output_ok$Pos+a
  output_ok$col<-b
  output_ok$col[output_ok$SNP %in% x$cof]<-'red'
  d<-(stats::aggregate(output_ok$xpos,list(output_ok$Chr),min)$x+stats::aggregate(output_ok$xpos,list(output_ok$Chr),max)$x)/2
  graphics::plot(output_ok$xpos,-log10(output_ok$pval),col=output_ok$col,pch=20,ylab=expression(-log[10](italic(p))),xaxt='n',xlab='chromosome',...)
  graphics::axis(1,tick=FALSE,at=d,labels=unique(output_ok$Chr))
  graphics::abline(h=x$bonf_thresh,lty=3,col='black')}

##' Plot
##'
##' Plot
##' @param x x
##' @param chrom chrom
##' @param pos1 pos1
##' @param pos2 pos2
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_region<-function(x,chrom,pos1,pos2){
  Chr<-Pos<-NULL # to avoid R CMD check issuing a NOTE
  region<-subset(x$out,Chr==chrom & Pos>=pos1 & Pos <=pos2)
  region$col<- if (chrom %% 2 == 0) {'gray60'} else {'gray10'}
  region$col[which(region$SNP %in% x$cof)]<-'red'
  graphics::plot(region$Pos,-log10(region$pval),type='p',pch=20,main=paste('chromosome',chrom,sep=''),xlab='position (bp)',ylab=expression(-log[10](italic(p))),col=region$col,xlim=c(pos1,pos2))
  graphics::abline(h=x$bonf_thresh,lty=3,col='black')}

##' Plot
##'
##' Plot
##' @param x x
##' @param step step
##' @param snp_info snp_info
##' @param pval_filt pval_filt
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, such as \code{main}
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_fwd_GWAS<-function(x,step,snp_info,pval_filt,...) {
  stopifnot(step<=length(x$pval_step))
  pval<-NULL # to avoid R CMD check issuing a NOTE
  output<-list(out=subset(merge(snp_info,x$pval_step[[step]]$out,by='SNP'),pval<=pval_filt),cof=x$pval_step[[step]]$cof,bonf_thresh=x$bonf_thresh)
  plot_GWAS(output,...)}

##' Plot
##'
##' Plot
##' @param x x
##' @param step step
##' @param snp_info snp_info
##' @param pval_filt pval_filt
##' @param chrom chrom
##' @param pos1 pos1
##' @param pos2 pos2
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_fwd_region<-function(x,step,snp_info,pval_filt,chrom,pos1,pos2) {
  stopifnot(step<=length(x$pval_step))
  pval<-NULL # to avoid R CMD check issuing a NOTE
  output<-list(out=subset(merge(snp_info,x$pval_step[[step]]$out,by='SNP'),pval<=pval_filt),cof=x$pval_step[[step]]$cof,bonf_thresh=x$bonf_thresh)
  plot_region(output,chrom,pos1,pos2)}

##' Plot
##'
##' Plot
##' @param x x
##' @param opt opt
##' @param snp_info snp_info
##' @param pval_filt pval_filt
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, such as \code{main}
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_opt_GWAS<-function(x,opt,snp_info,pval_filt,...) {
  pval<-NULL # to avoid R CMD check issuing a NOTE
  if (opt=='extBIC') {output<-list(out=subset(merge(snp_info,x$opt_extBIC$out,by='SNP'),pval<=pval_filt),cof=x$opt_extBIC$cof,bonf_thresh=x$bonf_thresh)
  plot_GWAS(output,...)}
  else if (opt=='mbonf') {output<-list(out=subset(merge(snp_info,x$opt_mbonf$out,by='SNP'),pval<=pval_filt),cof=x$opt_mbonf$cof,bonf_thresh=x$bonf_thresh)
  plot_GWAS(output,...)}
  else if (opt=='thresh') {output<-list(out=subset(merge(snp_info,x$opt_thresh$out,by='SNP'),pval<=pval_filt),cof=x$opt_thresh$cof,bonf_thresh=x$thresh)
  plot_GWAS(output,...)}
  else {cat('error! \n opt must be extBIC, mbonf or thresh')}}

##' Plot
##'
##' Plot
##' @param x x
##' @param opt opt
##' @param snp_info snp_info
##' @param pval_filt pval_filt
##' @param chrom chrom
##' @param pos1 pos1
##' @param pos2 pos2
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
plot_opt_region<-function(x,opt,snp_info,pval_filt,chrom,pos1,pos2) {
  pval<-NULL # to avoid R CMD check issuing a NOTE
  if (opt=='extBIC') {output<-list(out=subset(merge(snp_info,x$opt_extBIC$out,by='SNP'),pval<=pval_filt),cof=x$opt_extBIC$cof,bonf_thresh=x$bonf_thresh)
  plot_region(output,chrom,pos1,pos2)}
  else if (opt=='mbonf') {output<-list(out=subset(merge(snp_info,x$opt_mbonf$out,by='SNP'),pval<=pval_filt),cof=x$opt_mbonf$cof,bonf_thresh=x$bonf_thresh)
  plot_region(output,chrom,pos1,pos2)}
  else if (opt=='thresh') {output<-list(out=subset(merge(snp_info,x$opt_thresh$out,by='SNP'),pval<=pval_filt),cof=x$opt_thresh$cof,bonf_thresh=x$thresh)
  plot_region(output,chrom,pos1,pos2)}
  else {cat('error! \n opt must be extBIC, mbonf or thresh')}}

##' Plot
##'
##' Plot
##' @param x x
##' @param nsteps nsteps
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
qqplot_fwd_GWAS<-function(x,nsteps){
  stopifnot(nsteps<=length(x$pval_step))
  e<--log10(stats::ppoints(nrow(x$pval_step[[1]]$out)))
  ostep<-list()
  ostep[[1]]<--log10(sort(x$pval_step[[1]]$out$pval))
  for (i in 2:nsteps) {ostep[[i]]<--log10(sort(x$pval_step[[i]]$out$pval))}

  maxp<-ceiling(max(unlist(ostep)))

  graphics::plot(e,ostep[[1]],type='l',col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp))
  graphics::abline(0,1,col="dark grey")

  for (i in 2:nsteps) {
    graphics::par(new=T)
    graphics::plot(e,ostep[[i]],type='l',col=i,axes='F',xlab='',ylab='',xlim=c(0,max(e)+1),ylim=c(0,maxp))}
  graphics::legend(0,maxp,lty=1,pch=20,col=c(1:length(ostep)),paste(c(0:(length(ostep)-1)),'cof',sep=' '))
}

##' Plot
##'
##' Plot
##' @param x x
##' @param opt opt
##' @author V. Segura & B. J. Vilhjalmsson
##' @export
qqplot_opt_GWAS<-function(x,opt){
  if (opt=='extBIC') {
    e<--log10(stats::ppoints(nrow(x$opt_extBIC$out)))
    o<--log10(sort(x$opt_extBIC$out$pval))
    maxp<-ceiling(max(o))
    graphics::plot(e,o,type='l',col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),main=paste('optimal model according to extBIC'))
    graphics::abline(0,1,col="dark grey")}
  else if (opt=='mbonf') {
    e<--log10(stats::ppoints(nrow(x$opt_mbonf$out)))
    o<--log10(sort(x$opt_mbonf$out$pval))
    maxp<-ceiling(max(o))
    graphics::plot(e,o,type='l',col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),main=paste('optimal model according to mbonf'))
    graphics::abline(0,1,col="dark grey")}
  else if (opt=='thresh') {
    e<--log10(stats::ppoints(nrow(x$opt_thresh$out)))
    o<--log10(sort(x$opt_thresh$out$pval))
    maxp<-ceiling(max(o))
    graphics::plot(e,o,type='l',col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),main=paste('optimal model according to the user defined threshold'))
    graphics::abline(0,1,col="dark grey")}
  else {cat('error! \n opt must be extBIC, mbonf or thresh')}}
