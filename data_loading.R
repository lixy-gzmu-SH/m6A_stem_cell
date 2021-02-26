#basic loading
options(stringsAsFactors = F)
options(scipen = 200)
library(limma)
library(VennDiagram)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(edgeR)
library(ggrepel)
#library(ConsensusClusterPlus)
#library(survival)
#library(survminer)
library(reshape2)
#library(cluster)
library(ggsci)
library(Hmisc)
#library(ROCR)
#library(mice)
#library(pROC)
#library(glmnet)
#library(caret)
#library(papeR)
#library(ggthemes)
#library(risksetROC)
#library(rms)
library(clusterProfiler)
library(data.table)
#library(pathview)
#source("subtype_cluster.R")
source("ggcdf.R")
source("function.R")
source("pie_plot.R")
library(plyr)
library(ggrepel)
library(RColorBrewer)


# sampel
cell_line_sampe=c("HK1.E","HK1","N2.E","N2","N2.q","N2.q.E","NP460","NP460.E")
EBV_sample=cell_line_sampe[grep("E",cell_line_sampe)]
nEBV_sample=cell_line_sampe[grep("E",cell_line_sampe,invert = T)]

#color
tumor.color="#F7413C"
normal.color="#0005FC"
random.color="#808080"
tumor.color2="firebrick"
normal.color2="navy"
overlap.color="#A83690"

annotation_color1="#00B5BC"
annotation_color2="#F25B51"
annotation_color3="#F1EA0D"
annotation_color4="#058142"
annotation_color5="#F57F25"
annotation_color6="#D359A1"
annotation_color7="#3E6AB2"
annotation_color8="#9CBB5C"
annotation_color9="#8063A0"
annotation_color10="IndianRed1"
annotation_color11="Wheat"
annotation_color12="VioletRed"
annotation_color13="LightGreen"

#Peaks overlap
system("cat macs2/*.narrowPeak > macs2/peak_sample_peaks.bed")
system("intersectBed -a macs2/peak_sample_peaks.bed -b merge_peak/merged_Peak.bed -wa -wb > macs2/all_sample_in_merged_Peak.bed")
system("intersectBed -a macs2/all_sample_in_merged_Peak.bed -b merge_peak/all_metpeak.bed -wa -u > macs2/all_sample_in_macs2_and_metpeak.bed")

peak_in_macs2=system("awk -F \"\t\" '{print $14}' macs2/all_sample_in_merged_Peak.bed|sort|uniq|wc -l",intern = T)
#peak_in_metpeak=system("awk -F \"\t\" '{print $3}' merge_peak/all_metpeak.bed|sort|uniq|wc -l",intern = T)

peak_in_macs2=as.numeric(peak_in_macs2)
#peak_in_metpeak=as.numeric(peak_in_metpeak)

allpeak=fread("macs2/all_sample_in_macs2_and_metpeak.bed")
allpeak=data.frame(allpeak)
allpeak=cbind(allpeak,all=strsplit2(strsplit2(allpeak$V4,split="/")[,2],split="[.]")[,1])
allpeak$peak.id=paste(allpeak$V11,":",allpeak$V12,"-",allpeak$V13,sep="")
allpeak=allpeak[!duplicated(paste(allpeak$peak.id,allpeak$all,sep="_")),]

overlap_peak_bed=unique(allpeak$peak.id)
overlap_peak_bed=cbind(strsplit2(overlap_peak_bed,split=":|-"),overlap_peak_bed)
write.table(overlap_peak_bed,"annotation/overlap_peak.bed",sep="\t",quote = F,row.names = F,col.names = F)


#m6A peaks annotation
anno=read.table("annotation/human_EBV.anno.txt",header = F,sep="\t")
colnames(anno)=c("Chr","Start","End","Peak.id","Chr.gene","Start.gene","End.gene","Transcript.id",
                         "Nouse","Strand","Gene","Gene.type","Gene.site","Peak.position","Ensembl.gene.id","Level.1.gene.type",
                         "Level.2.gene.type")
rownames(anno)=anno$Peak.id

unanno=read.table("annotation/human_EBV.unanno.txt",header = F,sep="\t")
colnames(unanno)[4]="Peak.id"
unanno=merge(unanno,anno,by="Peak.id",all.x=T)
unanno=unanno[,colnames(anno)]
unanno[,c(1:3)]=strsplit2(unanno$Peak.id,split=":|-")
unanno[,c(5:7)]=unanno[,c(1:3)]

#EBV m6A peaks annotation by zhu kaiyu
EBV_annp=read.table("merge_peak/EBV_m6A_peaks_annotation.anno2.bed",sep="\t")

rownames(unanno)=unanno$Peak.id
unanno[paste(EBV_annp$V1,":",EBV_annp$V2,"-",EBV_annp$V3,sep=""),"Gene"]=gsub("gene_id ","",EBV_annp$V6)
unanno[paste(EBV_annp$V1,":",EBV_annp$V2,"-",EBV_annp$V3,sep=""),"Level.1.gene.type"]="EBV"
unanno[paste(EBV_annp$V1,":",EBV_annp$V2,"-",EBV_annp$V3,sep=""),"Level.2.gene.type"]="EBV"
unanno[paste(EBV_annp$V1,":",EBV_annp$V2,"-",EBV_annp$V3,sep=""),"Gene.site"]="EBV"

anno=rbind(anno,unanno)

anno[which(is.na(anno$Level.2.gene.type)),"Level.2.gene.type"]="Unknown"
rownames(anno)=anno$Peak.id

anno[which(anno$Gene.type=="coding" & anno$Level.2.gene.type!="mRNA"),"Level.2.gene.type"]="mRNA"
anno[which(anno$Gene.type=="coding" & anno$Level.2.gene.type!="mRNA"),"Level.1.gene.type"]="mRNA"
anno$Gene=gsub(" ","",anno$Gene)
rm(unanno)

#Find m6Am 5'UTR peaks
anno_5UTR=anno[which(anno$Gene.site=="5UTR"),]

gtf.temp=fread("/data/database/hg38/GENCODE/ebv-genome/human_EBV_Akata.gtf",sep="\t",skip = 5,data.table = F)
gtf.temp=cbind(gtf.temp,Transcript.id=strsplit2(strsplit2(gtf.temp$V9,split = "transcript_id ")[,2],split = ";")[,1])
gtf.temp$Transcript.id=gsub("\"","",gtf.temp$Transcript.id)
gtf.temp$Gene.site=gtf.temp$V3
gtf.temp$UTR.temp=strsplit2(gtf.temp[,9],split = ";")[,9]
gtf.temp[which(gtf.temp$Gene.site=="UTR" & gtf.temp$UTR.temp==" exon_number 1"),"Gene.site"]="5UTR"
gtf.temp[which(gtf.temp$Gene.site=="UTR" & gtf.temp$UTR.temp!=" exon_number 1"),"Gene.site"]="3UTR"
gtf.temp$UTR.temp=NULL

write.table(strsplit2(anno_5UTR$Peak.id,split=":|-"),"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")

#find m6Am motif
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -fo UTR5.peak.fa")
system("/data/software/homer2/bin/homer2 find -i UTR5.peak.fa -m BCA.motif -p 50 > BCA_peak_offset.txt")
BCA_in_5UTR_offset=read.table("BCA_peak_offset.txt",header=F)

anno_5UTR=overlap.anno[unique(BCA_in_5UTR_offset$V1),]
anno_5UTR=merge(gtf.temp,anno_5UTR,by="Transcript.id",all.y=T)
anno_5UTR=anno_5UTR[which(anno_5UTR$V3=="UTR"),]
anno_5UTR$Start=as.numeric(anno_5UTR$Start)
anno_5UTR$temp.start=anno_5UTR$V4-anno_5UTR$Start
anno_5UTR$temp.end=anno_5UTR$V5-anno_5UTR$Start
anno_5UTR[which(anno_5UTR$temp.start>0),"temp.start"]=1
anno_5UTR[which(anno_5UTR$temp.start<0),"temp.start"]=(-1)
anno_5UTR[which(anno_5UTR$temp.end>0),"temp.end"]=1
anno_5UTR[which(anno_5UTR$temp.end<0),"temp.end"]=(-1)
anno_5UTR=anno_5UTR[which((anno_5UTR$temp.start*anno_5UTR$temp.end)<=0),]
anno_5UTR.bed=cbind(anno_5UTR$V1,anno_5UTR$V4,anno_5UTR$V5,anno_5UTR$Peak.id,".",anno_5UTR$V7)
write.table(anno_5UTR.bed,"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -s -name -fo UTR5.peak.fa")

#get gene start with "A"
temp.utr5=read.table("UTR5.peak.fa",sep="\n")
temp.utr5=cbind(temp.utr5,substr(temp.utr5[,1],1,1))

i=2
n=nrow(temp.utr5)
temp.utr5=cbind(temp.utr5,type=NA)
while(i<=n){
  if(temp.utr5[i,2]=="A"){
    temp.utr5[c(i-1,i),"type"]="m6Am"
  }
  i=i+2
}
m6Am=na.omit(temp.utr5)
m6Am=m6Am[grep(">",m6Am$V1),]
m6Am=gsub(">","",m6Am$V1)
m6Am=strsplit2(m6Am,split="[(]")[,1]

anno=anno[setdiff(rownames(anno),m6Am),]

anno$Start=as.numeric(anno$Start)
anno$End=as.numeric(anno$End)


#all peak topology
all_topology=data.frame(cbind(anno[!is.na(anno$Gene),],type="All"))


#load gene expression matrix
gtf.gene=gtf.temp[which(gtf.temp$V3=="gene"),]
id2name=cbind(strsplit2(strsplit2(gtf.gene$V9,split="gene_id ")[,2],split = ";")[,1],
              strsplit2(strsplit2(gtf.gene$V9,split="gene_name ")[,2],split = ";")[,1])
id2name=gsub("\"","",id2name)
gtf_EBV=gtf.temp[grep("EBV",gtf.temp$V1),]
id2name_EBV=cbind(strsplit2(strsplit2(gtf_EBV$V9,split="gene_id ")[,2],split = ";")[,1])
id2name_EBV=gsub("\"","",id2name_EBV)
id2name=rbind(id2name,cbind(id2name_EBV,id2name_EBV))
rownames(id2name)=id2name[,1]

gene_count_data=read.table("gencode.rsem.count.txt",header=T,row.names = 1)
gene_count_data=cbind(gene=id2name[rownames(gene_count_data),2],gene_count_data)
gene_count_data[which(is.na(gene_count_data$gene)),"gene"]=rownames(gene_count_data[which(is.na(gene_count_data$gene)),])
gene_count_data=gene_count_data[!duplicated(gene_count_data$gene),]
rownames(gene_count_data)=gene_count_data$gene
gene_count_data$gene=NULL

gene_exp_data=read.table("gencode.rsem.TPM.txt",header=T,row.names = 1)
gene_exp_data=cbind(gene=id2name[rownames(gene_exp_data),2],gene_exp_data)
gene_exp_data[which(is.na(gene_exp_data$gene)),"gene"]=rownames(gene_exp_data[which(is.na(gene_exp_data$gene)),])
gene_exp_data=gene_exp_data[!duplicated(gene_exp_data$gene),]
rownames(gene_exp_data)=gene_exp_data$gene
gene_exp_data$gene=NULL

write.csv(gene_exp_data,"gene_exp_data.csv",quote = F,row.names = T,col.names = T)

#load peak enrichment matrix
rpkm_data=read.table("merge_peak/merge_Peak_abundance_normbytmn.txt",header=T,row.names = 1,fill=T)
rpkm_data=na.omit(rpkm_data)
rpkm_data=rpkm_data[anno$Peak.id,]
rpkm_ip_data=rpkm_data[,grep("ip",colnames(rpkm_data))]
rpkm_ip_data$HK1_3_2020.ip=rpkm_ip_data$HK1.ip
rpkm_input_data=rpkm_data[,grep("input",colnames(rpkm_data))]
rpkm_ip_norm=(rpkm_ip_data+1)/(rpkm_input_data+1)
colnames(rpkm_ip_norm)=gsub(".ip","",colnames(rpkm_ip_norm))
rpkm_ip_norm=na.omit(rpkm_ip_norm)

#WER
WER=c("METTL3","METTL14","METTL16","WTAP","RBM15","KIAA1429",
      "ALKBH5","FTO",
      "YTHDF1","YTHDF2","YTHDF3",
      "IGF2BP1","IGF2BP2","IGF2BP3",
      "YTHDC1","YTHDC2","SETD2","ZCCHC4","PRRC2A","DGCR8","CSTF2")



WE=c("METTL3","METTL14","METTL16","WTAP","RBM15","KIAA1429","ALKBH5","FTO","SETD2","CSTF2")
reader=c("YTHDF1","YTHDF2","YTHDF3","IGF2BP1","IGF2BP2","IGF2BP3","YTHDC1","YTHDC2")




RBP_POSTAR=read.table("all_RBP_from_POSTAR2")





