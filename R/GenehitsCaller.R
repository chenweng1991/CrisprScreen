#' GenehitsCaller
#'
#' This function allows you to collapse sgRNA into gene in order to identify gene hits.
#' @param Tier1Cut  The pvalue/FDR cutoff below which an sgRNA is considered tier1 sgRNA(very significant)
#' @param Tier2Cut  The pvalue/FDR cutoff between which to Tier1Cut an sgRNA is considered tier2 sgRNA(Somewhat significant)
#' @param Tier1Nmin Minimum number of supporting Tier1 sgRNA required to identify a gene hit (1,2,3,4,5,6)
#' @param TotalNmin Minimum number of supporting Tier1+Tier2 sgRNA required to identify a gene hit (1,2,3,4,5,6)
#' @param conflict  When getting rid of genes with conficting significant sgRNA. i,e. sgRNA show opposite direction. if "Onlytier1" then only consider confliction in Tier1 if "Alltiers" then consider both tier1 and tier2
#' @return this will return a list of two dataframe including gene hits that show the number of supporting sgRNA; and tier1+tier2 sgRNA data of all hits genes
#' @export
#' @examples
#' sgRNAdistribute(enrichment.data=P5.P6.enrichment)
GenehitsCaller<-function(enrichment.data=P5.P6.enrichment,Tier1Cut=0.001,Tier2Cut=0.01,Tier1Nmin=1,TotalNmin=2,conflict="Onlytier1"){
enrichment.data<-cbind(enrichment.data,Tiers=apply(enrichment.data,1,function(x){
if(as.numeric(as.character(x[5]))<=Tier1Cut){
  return("Tier1")
}else if(as.numeric(as.character(x[5]))<=Tier2Cut){
  return("Tier2")
}else{
return("")
}}))
enrichment.data.sig<-subset(enrichment.data,Tiers!="")
enrichment.data.sig$Gene<-as.character(enrichment.data.sig$Gene)
enrichment.data.sig$Tiers<-as.character(enrichment.data.sig$Tiers)
enrichment.sigGene<-table(enrichment.data.sig$Gene,enrichment.data.sig$Tiers)  %>% as.matrix(.)
enrichment.sigGene<-data.frame(Tier1=enrichment.sigGene[,1],Tier2=enrichment.sigGene[,2])
GeneHits<-subset(enrichment.sigGene,Tier1>=Tier1Nmin & (Tier1+Tier2)>=TotalNmin)
if(conflict=="Onlytier1"){
  subset(enrichment.data.sig,Gene %in% row.names(GeneHits)) %>% subset(.,Tiers=="Tier1") ->Hits.sig
}else if(conflict=="Alltiers"){
  subset(enrichment.data.sig,Gene %in% row.names(GeneHits)) -> Hits.sig
}
confictgenes<-c()
for (gene in unique(Hits.sig$Gene)){
  if(length(unique(subset(Hits.sig,Gene==gene)$Directions))==2){
    confictgenes<-c(confictgenes,gene)
  }
}
GeneHits<-GeneHits[!row.names(GeneHits) %in% confictgenes,]
GeneHits.sgRNA<-subset(Hits.sig,!Gene %in% confictgenes)
return(list(GeneHits=GeneHits,GeneHits.sgRNA=GeneHits.sgRNA))
}
