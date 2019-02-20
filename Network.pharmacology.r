#环境
#install.packages("readxl")
#install.packages('splitstackshape')
library(tidyr)
library(readxl)
library(splitstackshape)

#工作目录
setwd("C:/Users/summe/Desktop/Project/胆网络药理学/20190220 DATA")

#读取致病基因
gs<-read_excel("C0008350_disease_gda_summary (1).xlsx")
gs.name<-colnames(gs)
gs.gene<-gs[grep("symbol",gs.name)] 

#读取tcm有效kegg通路
tcm<-read_excel("batman-I2019-02-14-95328-1550137711-EnrichmentResult-Cutoff=20targets.xls")
tcm.all<-read.delim2("batman-I2019-02-14-95328-1550137711-cluster1-ScoreGT2targets.txt")

tcm.name<-colnames(tcm)
colnames(tcm)[6]<-"p-value|benjamini-value|coverage|EnrichRatio|gene names"
tcm.split<-strsplit(tcm$`p-value|benjamini-value|coverage|EnrichRatio|gene names`,split = "\\|")

tcm.gene<-c()

for (i in 1:length(tcm.split)) {
  tcm.gene.split<-strsplit(tcm.split[[i]][[5]],split = ";")
  tcm.gene<-append(tcm.gene,unlist(tcm.gene.split))
}

tcm.symbol<-data.frame("symbol"=tcm.gene)

#tcm调控和疾病共同基因

con.symbol<-merge(gs.gene,tcm.symbol)
con.symbol<-data.frame(table(con.symbol))
con.symbol<-con.symbol[1]
colnames(con.symbol)<-"symbol"
con.symbol.id<-merge(con.symbol,gs,by.x="symbol")

#go分析
#install.packages("colorspace")
#install.packages("stringi")
#source("http://bioconductor.org/biocLite.R")
#biocLite("DOSE")
#biocLite("clusterProfiler")
#biocLite("pathview")
library("clusterProfiler")
library("org.Hs.eg.db")

gs.go <- enrichGO(gene = con.symbol.id$geneid,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05)

#柱状图
barplot(gs.go, drop = TRUE, showCategory =20)

#点图
dotplot(gs.go,showCategory = 20)


#kegg分析
library("clusterProfiler")
gs.kegg<-enrichKEGG(gene = con.symbol.id$geneid, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)

#柱状图
barplot(gs.kegg, drop = TRUE, showCategory =20)

#点图
dotplot(gs.kegg,showCategory = 20)

#全部基因
tcm.all<-read.delim2("batman-I2019-02-14-95328-1550137711-cluster1-ScoreGT2targets.txt",header = FALSE)
compond.gene.all<-data.frame("compond"=NA,
                             "gene.name"=NA,
                             "gene.id"=NA)
for (i in 1:length(tcm.all)) {
  gene.local.name<-c()
  gene.local.id<-c()
  for (h in 3:dim(tcm.all)[2]){
  gene.s<-as.vector(tcm.all[i,h])
  gene.string.front<-regexpr("[|]",gene.s)
  gene.string.back<-regexpr("[(]",gene.s)
  gene.loci.id<-substr(gene.s,0,gene.string.front-1)
  gene.loci.name<-substr(gene.s,gene.string.front+1,gene.string.back-1)
  gene.local.name<-c(gene.local.name,gene.loci.name)
  gene.local.id<-c(gene.local.id,gene.loci.id)
}
  compond.gene<-data.frame( "compond"=as.vector(tcm.all[i,1]),
              "gene.name"=gene.local.name,
              "gene.id"=gene.local.id)
  compond.gene.all<-rbind(compond.gene.all,compond.gene)
}

#全部tcm调控和疾病共同基因
colnames(compond.gene.all)[2]<-"symbol"
compond.gene.all.symbol<-compond.gene.all[2]
con.symbol<-merge(gs.gene,compond.gene.all.symbol)

con.symbol<-data.frame(table(con.symbol))
con.symbol<-con.symbol[1]
colnames(con.symbol)<-"symbol"
con.symbol.id<-merge(con.symbol,gs,by.x="symbol")

#go分析
#install.packages("colorspace")
#install.packages("stringi")
#source("http://bioconductor.org/biocLite.R")
#biocLite("DOSE")
#biocLite("clusterProfiler")
#biocLite("pathview")
library("clusterProfiler")
library("org.Hs.eg.db")

gs.go <- enrichGO(gene = con.symbol.id$geneid,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05)

#柱状图
barplot(gs.go, drop = TRUE, showCategory =20)

#点图
dotplot(gs.go,showCategory = 20)


#kegg分析
library("clusterProfiler")
gs.kegg<-enrichKEGG(gene = con.symbol.id$geneid, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)

#柱状图
barplot(gs.kegg, drop = TRUE, showCategory =20)

#点图
dotplot(gs.kegg,showCategory = 20)
