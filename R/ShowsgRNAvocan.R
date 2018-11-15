#' ShowsgRNAvocano
#'
#' This function plot vocano for sgRNA and label the top sgRNAs
#' @param enrichment.data  The data prepared by "ComputeEnrichment"
#' @param pCut  p value/fdr cutoff for labeling sgRNAs
#' @param fcCut Foldchange cutoff for labeling sgRNAs
#' @param pdfname the file name for pdf of vocano plot
#' @param png,  if true then also make a png format. default is F
#' @param pngname the file name for png of vocano plot
#' @param pvalueadjusted To tell if the y axis used is the adjusted pvalue or original p value, default is TRUE
#' @return no return but generate a pdf or png
#' @export
#' @examples
#' ShowsgRNAvocano(enrichment.data=P5.P6.enrichment,pCut=1e-5,fcCut=2,pdfname="sgRNA.vocano.pdf",png=T,pngname="sgRNA.vocano.png")
ShowsgRNAvocano<-function(enrichment.data,pCut=1e-5,fcCut=2,pvalueadjusted=T,pdfname="sgRNA.vocano.pdf",png=F,pngname="sgRNA.vocano.png"){
require(ggrepel)
pdf(pdfname)
if(pvalueadjusted){
p<-ggplot(enrichment.data)+aes(log2(FCss),-log10(FDRs))+geom_point(size=0.5)+geom_text_repel(data=subset(enrichment.data,FDRs<pCut & abs(log2(FCss))>fcCut),aes(label=Gene),color="red")+theme(panel.background=element_blank(),axis.line=element_line(color="black",size=0.5))+ggtitle("sgRNA vocano plot")
}else{
p<-ggplot(enrichment.data)+aes(log2(FCss),-log10(pvalues))+geom_point(size=0.5)+geom_text_repel(data=subset(enrichment.data,FDRs<pCut & abs(log2(FCss))>fcCut),aes(label=Gene),color="red")+theme(panel.background=element_blank(),axis.line=element_line(color="black",size=0.5))+ggtitle("sgRNA vocano plot")
}
print(p)
dev.off()
if(png){
png(pngname)
print(p)
dev.off()
}
}
