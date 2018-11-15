
### Loading packages needed
```
library(dplyr)
library(qvalue)
library(ggplot2)
```

### Loading two readsnumber data for comparison
```
P5.readsnumber<-read.table("./data/F5_Ins-neg_pro-pos_ZF.readsnumber")
P6.readsnumber<-read.table("./data/F6_Ins-pos_pro-neg_ZF.readsnumber")
```
Readsnumber data looks like below

```
head(P5.readsnumber)
            V1           V2 V3
1          Pxk MGLibA_43941  6
2        Lrch1 MGLibB_29261 63
3 mmu-mir-8117 MGLibA_64983  8
4       Lrp2bp MGLibB_29338 97
5        Dyrk2 MGLibB_15213 22
6        S100b MGLibB_46842 17
```

### Merge the comparison data
- Here, we are going to merge two dataset together and add some dummy reads, default number is 5.
- Input two readsnumber data
- Give a sample name for each data

```
P5.P6.readsnumber<-PrepareData(readsnumber.1=P5.readsnumber,readsnumber.2=P6.readsnumber,name1="P5",name2="P6",dummyreads=5)
```

### Calculateing enrichment p value and fold change
- Here we are computing the p value and fold change of sgRNA enrichemnt in one of the two population
- In this version, a binomial test is applied.

```
P5.P6.enrichment<-ComputeEnrichment(data.merged=P5.P6.readsnumber,FDRmethod="qvalue")
```
### Quickly check out sgRNA p value and fold change distribution
```
sgRNAdistribute(enrichment.data=P5.P6.enrichment)
```

### Plot vocano plot for all sgRNA to show distribution and label the top sgRNA to have a sense if the experiment worked
```
ShowsgRNAvocano(enrichment.data=P5.P6.enrichment,pCut=1e-5,fcCut=1.5,pdfname="sgRNA.vocano.pdf",png=T,pngname="sgRNA.vocano.png")
```


### Collapse the sgRNA into gene level
In principle, we consider the p value/FDR and the nuber of supporting sgRNA and sgRNA direction confliction. To avoid overfitting, instead of calculating an overall score, we first classify sgRNAs into tiers. Tier1 is very strong sgRNA, tier2 is somewhat significant sgRNA. With the sgRNA classification, we set the number of sgRNA required in both Tier1 and tier2 to call a gene hit. TO control the quality of a gene hits list, we compute the globalfalse discovery rate by permutation

- The code below is to iteratively try different criteria combos to identify gene hits and compute gFDR. The purpose is to determine the most efficient parameter for gene hits calling
- If you already have the best parameter in mind, no need to run the code below

```
Criterias<-c()
gFDRs<-c()
i<-1
for (T1 in c(0.01,0.005,0.001,0.0001,0.00001)){  # to test the pvalue/FDR cutoff for tier1 sgRNA
  for(T2 in c(0.05,0.01,0.005)){  # to test the pvalue/FDR cutoff for tier2 sgRNA
    for(T1N in c(1,2)){  # to test hoiw many tier1 sgRNA required
      for(TTN in c(1,2,3,4)){ # to test how many tier1+tier2 sgRNA required
        for (CFL in c("Alltiers","Onlytier1")){ # to test to consider all tiers or only tier1 in ruling out conflicting genes
          print("###########")
          print(paste("#",i,"test"))
          i=i+1
          if (T2>T1){
            Cur.FDR<-ComputeGlobalFDR(enrichment.data=P5.P6.enrichment,Tier1Cut=T1,Tier2Cut=T2,Tier1Nmin=T1N,TotalNmin=TTN,conflict=CFL,shuffN=100)
            Criterias<-c(Criterias,paste(T1,T2,T1N,TTN,CFL,sep="_"))
            gFDRs<-rbind(gFDRs,Cur.FDR)
          }
        }
      }
    }
  }
}
row.names(gFDRs)<-Criterias
gFDRs<-data.frame(gFDRs,Criterias)
ggplot(gFDRs)+aes(hitGeneN.obs,gFDR)+geom_point()+geom_text_repel(data=subset(gFDRs,gFDR<0.1),aes(label=Criterias))  # TO visulize and find the best parameter combo
```

### Call gene hits using optimized parameter
Now with the best parameter in mind, we can call gene hits
- this will return a list of two dataframe including gene hits that show the number of supporting sgRNA; and tier1+tier2 sgRNA data of all hits genes
```
GenehitsCaller(enrichment.data=P5.P6.enrichment,Tier1Cut=0.001,Tier2Cut=0.01,Tier1Nmin=1,TotalNmin=2,conflict="Onlytier1")
```
