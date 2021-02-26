option(stringasFactor=FALSE)
options(stringsAsFactors = F)
library(data.table)
library(limma)
library(DESeq2)
library(edgeR)
library(pheatmap)
gene_count_data=read.table("gencode.rsem.count_v2.txt",
                           header = T,
                           row.names = 1)
gene2name=read.table("/data/xingyang/m6A_zhengjian/ensemble_id2gene_name.txt")
rownames(gene2name)=gene2name$V1
gene_count_data$Gene=gene2name[rownames(gene_count_data),"V2"]
gene_count_data=gene_count_data[!duplicated(gene_count_data$Gene),]
rownames(gene_count_data)=gene_count_data$Gene
gene_count_data$Gene=NULL

gene_exp_data=read.table("gencode.rsem.TPM_v2.txt",
                         header = T,
                         row.names = 1)
gene_exp_data$Gene=gene2name[rownames(gene_exp_data),"V2"]
gene_exp_data=gene_exp_data[!duplicated(gene_exp_data$Gene),]
rownames(gene_exp_data)=gene_exp_data$Gene
gene_exp_data$Gene=NULL

OM_vs_PM=diff.exp.gene(gene_count_data,
                        c(1:3),
                        c(4:6))

rownames(OM_vs_PM)=OM_vs_PM$Gene_ID

list_1=read.table("GO_REGULATION_OF_AUTOPHAGY")
list_2=read.table("GO_MITOPHAGY")

list_1_diff_table=OM_vs_PM[list_1$V1,]
list_2_diff_table=OM_vs_PM[list_2$V1,]

list_1_diff_table=na.omit(list_1_diff_table)
list_2_diff_table=na.omit(list_2_diff_table)


write.csv(list_1_diff_table,"GO_REGULATION_OF_AUTOPHAGY.csv",row.names = F,quote = F)
write.csv(list_2_diff_table,"GO_MITOPHAGY.csv",row.names = F,quote = F)
