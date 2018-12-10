#' sgRNAdistribute
#'
#' This function allows you to how the sgRNA distributed in terms of pvalue/fdr and fold change
#' @param enrichment.data  The data prepared by "ComputeEnrichment"
#' @param ptotest  a vector of p values to be tested
#' @param pnames  a vector of names of p values
#' @param FCtotest  a vector of fold change to be tested ( >1)
#' @param FCnames  a vector of names of fold change to be tested ( >1)
#' @param FC.lt  a vector of fold change to be tested ( <1)
#' @param FC.lt.names  a vector of names of fold change to be tested ( <1)
#' @return No return but print out useful information of p value and fdr distribution
#' @export
#' @examples
#' sgRNAdistribute(enrichment.data=P5.P6.enrichment,totest=c(0.05,0.01,0.001,0.0001,0.00001),pnames=c("<0.05","<0.01","<0.001","<0.0001","<0.00001"),FCtotest=c(1.2,1.5,2,5,10),FCnames=c(">1.2",">1.5",">2",">5",">10"),FC.lt=c(0.83,0.66,0.5,0.2,0.1),FC.lt.names=c("<0.83","<0.66","<0.5","<0.2","<0.1"))


sgRNAdistribute<-function(enrichment.data=P5.P6.enrichment,ptotest=c(0.05,0.01,0.001,0.0001,0.00001),pnames=c("<0.05","<0.01","<0.001","<0.0001","<0.00001"),FCtotest=c(1.2,1.5,2,5,10),FCnames=c(">1.2",">1.5",">2",">5",">10"),FC.lt=c(0.83,0.66,0.5,0.2,0.1),FC.lt.names=c("<0.83","<0.66","<0.5","<0.2","<0.1")){
sgRNAnumbers.p<-c()
for (pvalue in ptotest)
{
  sgRNAnumbers.p<-c(sgRNAnumbers.p,length(which(enrichment.data[,5]<=pvalue)))
}
names(sgRNAnumbers.p)<-pnames

sgRNAnumbers.F1<-c()
for (FC in FCtotest)
{
  sgRNAnumbers.F1<-c(sgRNAnumbers.F1,length(which(enrichment.data[,6]>FC)))
}
names(sgRNAnumbers.F1)<-FCnames

sgRNAnumbers.F2<-c()
for (FC in FC.lt)
{
  sgRNAnumbers.F2<-c(sgRNAnumbers.F2,length(which(enrichment.data[,6]<FC)))
}
names(sgRNAnumbers.F2)<-FC.lt.names
print(paste("totalo sgRNA number is",nrow(enrichment.data)))
print("sgRNA summarise by pvalue/fdr")
print(sgRNAnumbers.p)
print("sgRNA summarise by fold change")
print(sgRNAnumbers.F1)
print(sgRNAnumbers.F2)
}
