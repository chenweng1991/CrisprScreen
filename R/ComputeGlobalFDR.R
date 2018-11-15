#' ComputeGlobalFDR
#'
#' This function allows you compute the global FDR by permutation to evaluate the reliability of the gene hits. The parameter is exactly the same as GenehitsCaller
#' @param enrichment.data  The data prepared by "ComputeEnrichment"
#' @param Tier1Cut  The pvalue/FDR cutoff below which an sgRNA is considered tier1 sgRNA(very significant)
#' @param Tier2Cut  The pvalue/FDR cutoff between which to Tier1Cut an sgRNA is considered tier2 sgRNA(Somewhat significant)
#' @param Tier1Nmin Minimum number of supporting Tier1 sgRNA required to identify a gene hit (1,2,3,4,5,6)
#' @param TotalNmin Minimum number of supporting Tier1+Tier2 sgRNA required to identify a gene hit (1,2,3,4,5,6)
#' @param conflict  When getting rid of genes with conficting significant sgRNA. i,e. sgRNA show opposite direction. if "Onlytier1" then only consider confliction in Tier1 if "Alltiers" then consider both tier1 and tier2
#' @return This function will return a vector that include three component: expected gene hits number by permutation; observed gene hits nuber; global FDR
#' @export
#' @examples
#' sgRNAdistribute(enrichment.data=P5.P6.enrichment)
ComputeGlobalFDR<-function(enrichment.data=P5.P6.enrichment,Tier1Cut=0.001,Tier2Cut=0.01,Tier1Nmin=1,TotalNmin=2,conflict="Alltiers",shuffN=100){
Genehits.obs<-GenehitsCaller(enrichment.data=enrichment.data,Tier1Cut=Tier1Cut,Tier2Cut=Tier2Cut,Tier1Nmin=Tier1Nmin,TotalNmin=TotalNmin,conflict=conflict)
GenehitN.obs<-nrow(Genehits.obs$GeneHits)
GenehitN.shuffs<-c()
for (i in 1:shuffN){
print(i)
data.frame(Gene=sample(enrichment.data$Gene),enrichment.data[,2:ncol(enrichment.data)])  -> enrichment.data.shuf
GenehitN.shuff<-GenehitsCaller(enrichment.data=enrichment.data.shuf,Tier1Cut=Tier1Cut,Tier2Cut=Tier2Cut,Tier1Nmin=Tier1Nmin,TotalNmin=TotalNmin,conflict=conflict)$GeneHits %>% nrow
GenehitN.shuffs<-c(GenehitN.shuffs,GenehitN.shuff)
}
GenehitN.exp<-mean(GenehitN.shuffs)
return(c(hitGeneN.exp=GenehitN.exp,hitGeneN.obs=GenehitN.obs,gFDR=GenehitN.exp/GenehitN.obs))
}
