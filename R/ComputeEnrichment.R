
#' ComputeEnrichment
#'
#' This function allows you to compute enrichment fold change and p value(FDR) or each
#' @param data.merged  The data prepared by "PrepareData"
#' @param FDRmethod  The method used to adjust p value. "qvalue" or "bonferroni" or "none"
#' @return Return the data.merged that include p/fdr and fold change in the last column
#' @export
#' @examples
#' ComputeEnrichment(data.merged=P5.P6.readsnumber,FDRmethod="qvalue")
ComputeEnrichment<-function(data.merged=P5.P6.readsnumber,FDRmethod="qvalue"){
require(qvalue)
total.1<-sum(data.merged[,3][!is.na(data.merged[,3])])     #Total reads of third column, or readsnumber.1
total.2<-sum(data.merged[,4][!is.na(data.merged[,4])])     #Total reads of third column, or readsnumber.2
p.expected<-total.1/(total.1+total.2)    # expected ratio by total reads
pvalues<-apply(data.merged,1,function(x){
mod.cur<-binom.test(as.numeric(as.character(x[3])),(as.numeric(as.character(x[3]))+as.numeric(as.character(x[4]))),p.expected)
cur.pvalue<-mod.cur$p.value
return(cur.pvalue)
})
if(FDRmethod=="qvalue"){
FDR<-qvalue(pvalues)$qvalues
}else if(FDRmethod=="bonferroni"){
FDR<-p.adjust(pvalues,method="bonferroni")
}
if(FDRmethod!="none"){
data.merged<-cbind(data.merged,FDRs=FDR)
}else{
data.merged<-cbind(data.merged,pvalues=pvalues)
}
FCs<-apply(data.merged,1,function(x)
{
FC.cur<-(as.numeric(as.character(x[3]))/total.1)/(as.numeric(as.character(x[4]))/total.2)
return(FC.cur)
})
data.merged<-cbind(data.merged,FCss=FCs)
directions<-apply(data.merged,1,function(x){
if(as.numeric(as.character(x[6]))>1){
return("Negative")
}else{
return("Positive")
}})
data.merged<-cbind(data.merged,Directions=directions)
return(data.merged)
}
