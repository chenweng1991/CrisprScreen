
#' PrepareData
#'
#' This function allows you to prepare the data format for downstreme computation
#' @param readsnumber.1  The first readsnumber input file
#' @param readsnumber.2  The Second readsnumber input file
#' @param name1   Give a short name for data 1, such as "P5", "INShi", what ever distinguishable
#' @param name2   Give a short name for data 2, such as "P5", "INShi", what ever distinguishable
#' @param dummyreads Number of dummy reads to add on both in order to alleviate Low reads bias. Derfault is 5
#' @return return a formated merged data.frame with information from both experiment
#' @export
#' @examples
#' PrepareData(readsnumber.1=P5.readsnumber,readsnumber.2=P6.readsnumber,name1="P5",name2="P6",dummyreads=5)
PrepareData<-function(readsnumber.1=P5.readsnumber,readsnumber.2=P6.readsnumber,name1="P5",name2="P6",dummyreads=5)
{
data.merged<-merge(readsnumber.1,readsnumber.2,by="V2",all=T)
gene<-apply(data.merged,1,function(x){if(is.na(x[2]))return(x[4])else return(x[2]) })
data.merged<-cbind(gene,data.merged[,c(1,3,5)])
colnames(data.merged)<-c("Gene","sgRNAid",name1,name2)
data.merged[is.na(data.merged)]<-0
cbind(data.merged[,1:2],(data.merged[,3:4]+dummyreads)) ->data.merged
return(data.merged)
}
