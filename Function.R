#drop duplicated peaks
drop_duplicated_peak=function(temp){
  temp=temp[!duplicated(temp[,14]),]
  return(temp)
}



#specific peaks associate with clinic factor
peakAclinic=function(peakandfreq,clinic_data){
  peakandfreq=data.frame(peakandfreq)
  i=1
  clinic_factor=c("gender","age","fenhua_class","smoke","drink","nerve",
                  "xueguan","node","TNM_3")
  unique_peak=unique(peakandfreq$Peak.id)
  n=length(clinic_factor)
  result=c()
  while(i<=n){
    factor.type=unique(clinic_data[,clinic_factor[i]])
    factor.type=factor.type[order(factor.type,decreasing = T)]
    factor.name=clinic_factor[i]
    j=1
    k=length(factor.type)
    print(factor.type)
    clinic_factor_sample=c()
    while(j<=k){
      sample2group=cbind(Sample=clinic_data[which(clinic_data[,clinic_factor[i]]==factor.type[j]),"t_id"],
                         Group=factor.type[j])
      clinic_factor_sample=rbind(clinic_factor_sample,sample2group)
      j=j+1
    }
    j=1
    k=length(unique_peak)
    while(j<=k){
      support_sample=unique(peakandfreq[which(peakandfreq$Peak.id==unique_peak[j]),'Sample'])
      if(length(factor.type)>2){
        p=1
        q=length(factor.type)
        while(p<=q){
          peak_freq_matrix=matrix(
            c(table(merge(data.frame(Sample=support_sample),
                          clinic_factor_sample,by="Sample")[,"Group"])[as.character(factor.type[p])],
              table(merge(data.frame(Sample=setdiff(gsub("X","",tumor_sample_in_peak),support_sample)),
                          clinic_factor_sample,by="Sample")[,"Group"])[as.character(factor.type[p])],
              sum(na.omit(as.numeric(table(merge(data.frame(Sample=support_sample),
                                                 clinic_factor_sample,by="Sample")[,"Group"])[as.character(factor.type[-p])]))),
              sum(na.omit(as.numeric(table(merge(data.frame(Sample=setdiff(gsub("X","",tumor_sample_in_peak),support_sample)),
                                                 clinic_factor_sample,by="Sample")[,"Group"])[as.character(factor.type[-p])])))),2,2)
          peak_freq_matrix[which(is.na(peak_freq_matrix))]=0
          pvalue.greater=fisher.test(peak_freq_matrix,alternative = "greater")
          pvalue.greater=pvalue.greater$p.value
          pvalue.less=fisher.test(peak_freq_matrix,alternative = "less")
          pvalue.less=pvalue.less$p.value
          difftype=c("Greater","Less")[which.min(c(pvalue.greater,pvalue.less))]
          result=rbind(result,
                       c(unique_peak[j],paste(clinic_factor[i],": ",factor.type[p]," vs Rest",sep=""),difftype,
                         c(pvalue.greater,pvalue.less)[which.min(c(pvalue.greater,pvalue.less))],as.numeric(peak_freq_matrix)))
          p=p+1
        }
        
      } else {
        peak_freq_matrix=matrix(
          c(table(merge(data.frame(Sample=support_sample),clinic_factor_sample,by="Sample")[,"Group"])[as.character(factor.type)],
            table(merge(data.frame(Sample=setdiff(gsub("X","",tumor_sample_in_peak),support_sample)),clinic_factor_sample,by="Sample")[,"Group"])[as.character(factor.type)]),2,2)
        peak_freq_matrix[which(is.na(peak_freq_matrix))]=0
        pvalue.greater=fisher.test(peak_freq_matrix,alternative = "greater")
        pvalue.greater=pvalue.greater$p.value
        pvalue.less=fisher.test(peak_freq_matrix,alternative = "less")
        pvalue.less=pvalue.less$p.value
        difftype=c("Greater","Less")[which.min(c(pvalue.greater,pvalue.less))]
        result=rbind(result,
                     c(unique_peak[j],clinic_factor[i],difftype,
                       c(pvalue.greater,pvalue.less)[which.min(c(pvalue.greater,pvalue.less))],as.numeric(peak_freq_matrix)))
      }
      
      j=j+1
    }
    i=i+1
  }
  result=data.frame(result)
  colnames(result)[1:4]=c("Peak.id","Clinic.factor","Difftype","P.value")
  return(result)
}

#topology plot
ggtopo=function(anno_table){
  anno_table=anno_table[which(anno_table$Gene.type=="coding"),]
  anno_table=anno_table[which(anno_table$Gene.site!="intron"),]
  anno_table[which(anno_table$Gene.site=="CDS"),'Peak.position']=anno_table[which(anno_table$Gene.site=="CDS"),'Peak.position']+100
  anno_table[which(anno_table$Gene.site=="3UTR"),'Peak.position']=anno_table[which(anno_table$Gene.site=="3UTR"),'Peak.position']+200
  ggplot(anno_table)+geom_line(aes(x=anno_table$Peak.position,color=anno_table$type),stat = "density",size=1.5)+
    geom_vline(aes(xintercept=100), colour="#40004B", linetype="dashed")+
    geom_vline(aes(xintercept=200), colour="#40004B", linetype="dashed")+
    scale_x_continuous(breaks=c(50, 150, 250), labels=c("5UTR", "CDS", "3UTR"),limits = c(0,300))+
    xlab("")+ylab("Density")+
    theme_classic()+theme(axis.text.x = element_text(size=12),legend.title = element_blank(),
                          legend.text = element_text(size = 7),legend.justification=c(1,1), legend.position=c(0.3,0.99))
}

ggtopo_for_count=function(anno_table){
  anno_table=anno_table[which(anno_table$Gene.type=="coding"),]
  anno_table=anno_table[which(anno_table$Gene.site!="intron"),]
  anno_table[which(anno_table$Gene.site=="CDS"),'Peak.position']=anno_table[which(anno_table$Gene.site=="CDS"),'Peak.position']+100
  anno_table[which(anno_table$Gene.site=="3UTR"),'Peak.position']=anno_table[which(anno_table$Gene.site=="3UTR"),'Peak.position']+200
  ggplot(anno_table)+geom_line(aes(x=anno_table$Peak.position,color=anno_table$type),stat = "count",size=1.2)+
    geom_vline(aes(xintercept=100), colour="#40004B", linetype="dashed")+
    geom_vline(aes(xintercept=200), colour="#40004B", linetype="dashed")+
    scale_x_continuous(breaks=c(50, 150, 250), labels=c("5UTR", "CDS", "3UTR"))+
    xlab("")+ylab("Count")+
    theme_classic()+theme(axis.text.x = element_text(size=12),legend.title = element_blank(),
                          legend.text = element_text(size = 7),legend.justification=c(1,1), legend.position=c(0.99,0.9))
}




#clinic_heatmap
clinic_heatmap=function(peak_associate_with_clinic,top=10,peak2gene){
  peak_associate_with_clinic$P.value=as.numeric(peak_associate_with_clinic$P.value)
  peak_associate_with_clinic$log2FC=as.numeric(peak_associate_with_clinic$log2FC)
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$P.value<0.05),]
  peak_associate_with_clinic=merge(peak_associate_with_clinic,peak2gene,by="Peak.id",all.x = T)
  peak_associate_with_clinic=peak_associate_with_clinic[which(!is.na(peak_associate_with_clinic$gene)),]
  peak_associate_with_clinic=as.data.frame(cbind(peak_associate_with_clinic,type="non"))
  peak_associate_with_clinic=peak_associate_with_clinic[order(peak_associate_with_clinic$P.value),]
  i=1
  Clinic.factor=unique(peak_associate_with_clinic$Clinic.factor)
  n=length(Clinic.factor)
  while(i<=n){
    if(length(peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'Peak.id'])<top){
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'type']='sig'
    } else {
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i])[1:top],'type']='sig'
    }
    i=i+1
  }
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$type=="sig"),]
  peak_associate_with_clinic$P.value=as.numeric(peak_associate_with_clinic$P.value)
  peak_associate_with_clinic$Peak.id=paste(peak_associate_with_clinic$Peak.id,":",peak_associate_with_clinic$gene,sep="")
  peak_associate_with_clinic$Peak.id=ordered(peak_associate_with_clinic$Peak.id,
                                             levels=names(table(peak_associate_with_clinic$Peak.id)[
                                               order(table(peak_associate_with_clinic$Peak.id),decreasing = T)]))
  p=ggplot(peak_associate_with_clinic)+
    geom_point(aes(y=peak_associate_with_clinic$Clinic.factor,
                   x=peak_associate_with_clinic$Peak.id,size=(-log10(peak_associate_with_clinic$P.value)),
                   color=as.factor(peak_associate_with_clinic$Difftype)))+
    scale_color_manual(values=c("#d50709","#5707d5"),name="log2(Flod change)")+
    scale_size_continuous(name="-log10(P value)")+
    xlab("Clinic factor")+ylab("Tumor specific peaks top10 in each clinic factor")+
    theme_classic()+theme(axis.text.x = element_text(angle=90))
  return(p)
}

clinic_heatmap_for_common_hyper=function(peak_associate_with_clinic,top=10,peak2gene){
  peak_associate_with_clinic$FDR=as.numeric(peak_associate_with_clinic$FDR)
  peak_associate_with_clinic$log2FC=as.numeric(peak_associate_with_clinic$log2FC)
  peak_associate_with_clinic$P.value=as.numeric(peak_associate_with_clinic$P.value)
  #peak_associate_with_clinic$Contribution=as.numeric(peak_associate_with_clinic$Contribution)
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$FDR<0.1),]
  peak_associate_with_clinic=merge(peak_associate_with_clinic,peak2gene,by="Peak.id",all.x = T)
  peak_associate_with_clinic=peak_associate_with_clinic[which(!is.na(peak_associate_with_clinic$gene)),]
  peak_associate_with_clinic=as.data.frame(cbind(peak_associate_with_clinic,type="non"))
  peak_associate_with_clinic=peak_associate_with_clinic[order(abs(peak_associate_with_clinic$log2FC),decreasing = T),]
  peak_associate_with_clinic=peak_associate_with_clinic[order(peak_associate_with_clinic$FDR),]
  
  i=1
  Clinic.factor=unique(peak_associate_with_clinic$Clinic.factor)
  n=length(Clinic.factor)
  while(i<=n){
    if(nrow(peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),])<1){
      i=i+1
      next
    }
    if(length(peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'Peak.id'])<=top){
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'type']='sig'
    } else {
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i])[1:top],'type']='sig'
    }
    i=i+1
  }
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$type=="sig"),]
  peak_associate_with_clinic$P.value=as.numeric(peak_associate_with_clinic$P.value)
  peak_associate_with_clinic$Peak.id=paste(peak_associate_with_clinic$Peak.id,":",peak_associate_with_clinic$gene,sep="")
  peak_associate_with_clinic$Peak.id=ordered(peak_associate_with_clinic$Peak.id,
                                             levels=names(table(peak_associate_with_clinic$Peak.id))[
                                               order(table(peak_associate_with_clinic$Peak.id),decreasing = T)])
  peak_associate_with_clinic=data.frame(cbind(peak_associate_with_clinic,Difftype="Up"))
  peak_associate_with_clinic[which(peak_associate_with_clinic$log2FC<0),"Difftype"]="Down"
  #peak_associate_with_clinic$Clinic.factor=ordered(peak_associate_with_clinic$Clinic.factor,
  #                                                 levels=rev(c("age","gender","drink","smoke","TNM_3","xueguan","node","nerve",
  #                                                          "fenhua_class")))
  peak_associate_with_clinic$Clinic.factor=ordered(peak_associate_with_clinic$Clinic.factor,
                                                   levels=rev(c("age","gender","drink","smoke","TNM_3","xueguan","node","nerve",
                                                                "fenhua_class : 1  vs Rest","fenhua_class : 2  vs Rest",
                                                                "fenhua_class : 3  vs Rest")))
  
  write.table(sort(peak_associate_with_clinic$Peak.id),"analysis/new_figure/peak_ass.temp.txt",row.names = F,
              quote = F)
  p=ggplot(peak_associate_with_clinic)+
    geom_point(aes(y=peak_associate_with_clinic$Clinic.factor,
                   x=peak_associate_with_clinic$Peak.id,size=((abs(peak_associate_with_clinic$log2FC))),
                   color=as.factor(peak_associate_with_clinic$Difftype)))+
    scale_color_manual(values=c(Down="#5707d5",Up="#d50709"),name="log2FC")+
    scale_size_continuous(name="log2(FC)")+
    xlab("Peak")+ylab("Top10 in each clinic factor")+
    theme_classic()+theme(axis.text.x = element_text(angle=90))
  return(p)
}

#shapiro test for row
row.shapiro.test=function(common_peak){
  i=1
  n=nrow(common_peak)
  test=c()
  while(i<=n){
    pvalue=shapiro.test(as.numeric(common_peak[i,]))
    test=rbind(test,c(rownames(common_peak)[i],pvalue$p.value,if(pvalue$p.value>0.1){"yes"}else{"no"}))
    i=i+1
  }
  test=data.frame(test)
  colnames(test)=c("Peak.id","P.value","norm.type")
  return(test)
  
}

#wilcox test for tumor and normal
row.wilcox.test=function(common_peak,tumor_sample_in_peak,normal_sample_in_peak,paired.set=F){
  i=1
  n=nrow(common_peak)
  test=c()
  while(i<=n){
    pvalue.paired=wilcox.test(as.numeric(common_peak[i,tumor_sample_in_peak]),
                              as.numeric(common_peak[i,normal_sample_in_peak]),paired = paired.set)
    log2FC.paired=log2(mean(as.numeric(common_peak[i,tumor_sample_in_peak]))/
                         mean(as.numeric(common_peak[i,normal_sample_in_peak])))
    test=rbind(test,c(rownames(common_peak)[i],pvalue.paired$p.value,log2FC.paired))
    i=i+1
  }
  test=data.frame(test)
  colnames(test)=c("Peak.id","P.value.paired","log2FC.paired")
  test=cbind(test,FDR.paired=p.adjust(test$P.value.paired,method = "fdr"))
  test$P.value.paired=as.numeric(test$P.value.paired)
  test$log2FC.paired=as.numeric(test$log2FC.paired)
  return(test)
  
}


#different expression
diff.exp.gene=function(countmatrix,group_1,group_2){
  pheno=rbind(cbind(Status=0,sample=colnames(countmatrix)[group_1]),
              cbind(Status=1,sample=colnames(countmatrix)[group_2]))
  
  countmatrix=as.matrix(countmatrix)
  countmatrix=round(countmatrix)
  
  #                                                             DESeq2
  cond <- factor(pheno[,1])
  dds <- DESeqDataSetFromMatrix(countmatrix, DataFrame(cond), ~ cond)
  dds <- DESeq(dds)
  res <-  results(dds,independentFiltering = F)
  resOrdered <- res[order(res$log2FoldChange,decreasing = T),]
  resOrdered=as.data.frame(resOrdered)
  resOrdered=data.frame("Gene_ID"=rownames(resOrdered),resOrdered)
  
  degene_deseq=resOrdered
  print(dim(degene_deseq))
  #                                                             edgeR
  d <- DGEList(counts=countmatrix, group=cond)
  d <- calcNormFactors(d) 
  design <- model.matrix(~cond)
  d <- estimateDisp(d,design)
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit, coef=2) 
  print(dim(lrt$table))
  lrt=data.frame(Gene_ID=rownames(lrt$table),lrt$table)
  lrt=data.frame(cbind(lrt,FDR=p.adjust(lrt$PValue,method = "fdr")))
  degene_edger=lrt
  print(dim(degene_edger))
  co_diff=merge(degene_deseq,degene_edger,all=T,by="Gene_ID")
  #co_diff=na.omit(co_diff)
  return(co_diff)
}


#pathway barplot
pathway_barplot_for_GO=function(GO,color="#7292C5"){
  GO$PValue=as.numeric(GO$PValue)
  logp=-log10(GO$PValue)
  GO$Term=as.character(GO$Term)
  #temp=strsplit(GO$Term,"~")
  #temp=as.data.frame(temp)
  #temp=t(temp)
  #bar_plot=cbind(temp[,2],logp)
  bar_plot=cbind(GO$Term,logp)
  bar_plot=as.data.frame(bar_plot)
  bar_plot[,2]=as.numeric(as.character(bar_plot[,2]))
  p=ggplot(bar_plot,aes(reorder(x=bar_plot[,1],logp),y=bar_plot[,2]))+
    geom_bar(stat="identity",fill=color,width = 0.5)+coord_flip()+
    xlab("")+ylab("GO pathway enrichment -log10(p value)")+
    theme(axis.text.y = element_text(size=12))+
    theme_classic()
  return(p)
}

pathway_barplot_for_KEGG=function(GO,color="#7292C5"){
  GO$PValue=as.numeric(GO$PValue)
  logp=-log10(GO$PValue)
  GO$Term=as.character(GO$Term)
  #temp=strsplit(GO$Term,":")
  #temp=as.data.frame(temp)
  #temp=t(temp)
  #bar_plot=cbind(temp[,2],logp)
  bar_plot=cbind(GO$Term,logp)
  bar_plot=as.data.frame(bar_plot)
  bar_plot[,2]=as.numeric(as.character(bar_plot[,2]))
  p=ggplot(bar_plot,aes(reorder(x=bar_plot[,1],logp),y=bar_plot[,2]))+
    geom_bar(stat="identity",fill=color,width = 0.5)+coord_flip()+
    xlab("")+ylab("KEGG pathway enrichment -log10(p value)")+
    theme(axis.text.y = element_text(size=12))+theme_classic()
  return(p)
}

#common hyper associate with clinic
peak_median_for_clinic=function(peak_for_analysis,peak_matrix,clinic_data){
  i=1
  clinic_factor=c("gender","age","fenhua_class","smoke","drink","nerve",
                  "xueguan","node","TNM_3")
  n=length(clinic_factor)
  result=c()
  while(i<=n){
    factor.type=unique(clinic_data[,clinic_factor[i]])
    factor.type=factor.type[order(factor.type,decreasing = T)]
    factor.name=clinic_factor[i]
    j=1
    k=length(factor.type)
    print(factor.type)
    clinic_factor_sample=c()
    while(j<=k){
      sample2group=cbind(Sample=clinic_data[which(clinic_data[,clinic_factor[i]]==factor.type[j]),"t_id"],Group=factor.type[j])
      clinic_factor_sample=rbind(clinic_factor_sample,sample2group)
      j=j+1
    }
    clinic_factor_sample=data.frame(clinic_factor_sample)
    j=1
    k=length(peak_for_analysis)
    while(j<=k){
      p=1
      q=length(factor.type)
      if(q>2){
        while(p<=q){
          group_p=as.numeric(peak_matrix[peak_for_analysis[j],
                                         intersect(tumor_sample_in_peak,
                                                   paste(sep="","X",clinic_factor_sample[
                                                     which(clinic_factor_sample$Group==factor.type[p]),"Sample"]))])
          group_non_p=as.numeric(peak_matrix[peak_for_analysis[j],
                                             intersect(tumor_sample_in_peak,
                                                       paste(sep="","X",clinic_factor_sample[
                                                         which(clinic_factor_sample$Group!=factor.type[p]),"Sample"]))])
          pvalue=wilcox.test(group_p,group_non_p)
          pvalue=pvalue$p.value
          log2FC=log2(mean(group_p)/mean(group_non_p))
          result=rbind(result,
                       c(peak_for_analysis[j],paste(clinic_factor[i],":",factor.type[p]," vs Rest"),pvalue,log2FC))
          p=p+1
        } 
      } else {
        group_p=as.numeric(peak_matrix[peak_for_analysis[j],
                                       intersect(tumor_sample_in_peak,
                                                 paste(sep="","X",clinic_factor_sample[
                                                   which(clinic_factor_sample$Group==factor.type[1]),"Sample"]))])
        group_non_p=as.numeric(peak_matrix[peak_for_analysis[j],
                                           intersect(tumor_sample_in_peak,
                                                     paste(sep="","X",clinic_factor_sample[
                                                       which(clinic_factor_sample$Group!=factor.type[1]),"Sample"]))])
        pvalue=wilcox.test(group_p,group_non_p)
        pvalue=pvalue$p.value
        log2FC=log2(mean(group_p)/mean(group_non_p))
        result=rbind(result,
                     c(peak_for_analysis[j],clinic_factor[i],pvalue,log2FC))
        
      }
      j=j+1
    }
    i=i+1
  }
  result=data.frame(result)
  colnames(result)=c("Peak.id","Clinic.factor","P.value","log2FC")
  return(result)
}



#tumor specific peak survival
specific_peak_survival=function(peak_for_analysis,peak_freq,clinic_data_for_analysis,tumor_sample_in_peak){
  i=1
  n=length(peak_for_analysis)
  rownames(clinic_data_for_analysis)=clinic_data_for_analysis$t_id
  clinic_data_for_analysis=clinic_data_for_analysis[gsub("X","",tumor_sample_in_peak),]
  surv=Surv(time = clinic_data_for_analysis$month,clinic_data_for_analysis$status==1)
  result=c()
  while(i<=n){
    peak.sample=data.frame(t_id=unique(peak_freq[which(peak_freq$Peak.id==peak_for_analysis[i]),"Sample"]),type="Have peaks")
    print(peak.sample)
    peak_factor=data.frame(t_id=clinic_data_for_analysis$t_id)
    peak_factor=merge(peak_factor,peak.sample,by="t_id",all.x=T)
    peak_factor[which(is.na(peak_factor$type)),"type"]="Have not peaks"
    #peak_factor=peak_factor[!duplicated(peak_factor$t_id),]
    rownames(peak_factor)=peak_factor$t_id
    peak_factor[clinic_data_for_analysis$t_id,]
    fit=coxph(surv~peak_factor$type)
    fit=summary(fit)
    pvalue=fit$sctest[3]
    HR=fit$conf.int[1]
    HR.lower95=fit$conf.int[3]
    HR.upper95=fit$conf.int[4]
    result=rbind(result,c(peak_for_analysis[i],pvalue,HR,HR.lower95,HR.upper95))
    i=i+1
  }
  result=data.frame(result)
  colnames(result)=c("Peak.id","P.Value","HR","HR.lower95","HR.upper95")
  return(result)
}

specific_peak_survival_for_common_hyper=function(peak_for_analysis,clinic_data_for_analysis,tumor_sample_in_peak,peak_matrix){
  i=1
  n=length(peak_for_analysis)
  rownames(clinic_data_for_analysis)=clinic_data_for_analysis$t_id
  clinic_data_for_analysis=clinic_data_for_analysis[gsub("X","",tumor_sample_in_peak),]
  surv=Surv(time = clinic_data_for_analysis$month,clinic_data_for_analysis$status==1)
  result=c()
  while(i<=n){
    pvalue=c()
    peak.value=peak_matrix[peak_for_analysis[i],]
    if(median(as.numeric(peak.value))==0){
      peak.lowMethylation=peak.value[which(peak.value==median(as.numeric(peak.value)))]
      peak.highMethylation=peak.value[which(peak.value>median(as.numeric(peak.value)))]
    }else {
      peak.lowMethylation=peak.value[which(peak.value<median(as.numeric(peak.value)))]
      peak.highMethylation=peak.value[which(peak.value>=median(as.numeric(peak.value)))]
      peak_factor=rbind(cbind(t_id=gsub("X","",names(peak.highMethylation)),type="High methylation"),
                        cbind(t_id=gsub("X","",names(peak.lowMethylation)),type="Low methylation"))
    }
    peak_factor=rbind(cbind(t_id=gsub("X","",names(peak.highMethylation)),type="High methylation"),
                      cbind(t_id=gsub("X","",names(peak.lowMethylation)),type="Low methylation"))
    
    peak_factor=data.frame(peak_factor)
    rownames(peak_factor)=peak_factor$t_id
    peak_factor=peak_factor[clinic_data_for_analysis$t_id,]
    peak_factor$type=ordered(peak_factor$type,levels=c("Low methylation","High methylation"))
    fit=coxph(surv~peak_factor$type)
    fit=summary(fit)
    pvalue=fit$sctest[3]
    HR=fit$conf.int[1]
    HR.lower95=fit$conf.int[3]
    HR.upper95=fit$conf.int[4]
    HR.pvlue=fit$coefficients[5]
    result=rbind(result,c(peak_for_analysis[i],pvalue,HR,HR.lower95,HR.upper95,HR.pvlue))
    i=i+1
  }
  result=data.frame(result)
  colnames(result)=c("Peak.id","LogRank.P","HR","HR.lower95","HR.upper95","HR.Pvalue")
  result$LogRank.P=as.numeric(result$LogRank.P)
  result$HR=as.numeric(result$HR)
  result$HR.lower95=as.numeric(result$HR.lower95)
  result$HR.upper95=as.numeric(result$HR.upper95)
  result$HR.Pvalue=as.numeric(result$HR.Pvalue)
  return(result)
}


#survival for plot
survival_plot=function(peak_for_analysis,clinic_data_for_analysis,tumor_sample_in_peak,peak_matrix,label,output_filename,anno){
  
}




#heatmap for subtype
heatmap_for_subtype=function(subtype_sample_id,peak_matrix,TCGA_subtype,clinic_data_for_analysis,output_dir,draw_peak){
  peak_matrix=peak_matrix[draw_peak,]
  jj=1
  number_cluster=unique(subtype_sample_id$subtype)
  depeak_subtype=c()
  while(jj<=length(number_cluster)){
    group_p=rownames(subtype_sample_id[which(subtype_sample_id$subtype==unique(subtype_sample_id$subtype)[jj]),])
    group_non_p=rownames(subtype_sample_id[which(subtype_sample_id$subtype!=unique(subtype_sample_id$subtype)[jj]),])
    #subtype heatmap
    i=1
    n=nrow(peak_matrix)
    depeak_temp=c()
    while(i<=n){
      pvalue=wilcox.test(as.numeric(peak_matrix[i,group_p]),as.numeric(peak_matrix[i,group_non_p]))
      pvalue=pvalue$p.value
      log2FC=log2(mean(as.numeric(peak_matrix[i,group_p]))/mean(as.numeric(peak_matrix[i,group_non_p])))
      depeak_temp=rbind(depeak_temp,c(rownames(peak_matrix)[i],log2FC,pvalue))
      i=i+1
    }
    depeak_temp=data.frame(depeak_temp)
    depeak_temp$X2=as.numeric(depeak_temp$X2)
    depeak_temp$X3=as.numeric(depeak_temp$X3)
    depeak_temp=data.frame(cbind(depeak_temp,FDR=p.adjust(depeak_temp$X3,method="fdr")))
    if(unique(subtype_sample_id$subtype)[jj]=="s1"){
      depeak_temp=depeak_temp[which(depeak_temp$FDR<0.01 & (depeak_temp$X2)<(-0.58)),]
    }else{
      depeak_temp=depeak_temp[which(depeak_temp$FDR<0.01 & (depeak_temp$X2)>0.58),]
      
    }
    if(nrow(depeak_temp)<1){
      return(print("No peak in a subtype!"))
    }
    depeak_temp=depeak_temp[order(depeak_temp$X2,decreasing = T),]
    depeak_subtype=rbind(depeak_subtype,cbind(depeak_temp,Subtype=unique(subtype_sample_id$subtype)[jj]))
    print(jj)
    jj=jj+1
  }
  depeak_subtype=data.frame(depeak_subtype)
  colnames(depeak_subtype)=c("Peak.id","log2FC","P.value","FDR","Subtype")
  depeak_subtype=depeak_subtype[order(abs(depeak_subtype$log2FC),decreasing = T),]
  depeak_subtype=depeak_subtype[!duplicated(depeak_subtype$Peak.id),]
  depeak_subtype=depeak_subtype[order(depeak_subtype$Subtype),]
  print("order the data")
  de_peak_exp_plot=peak_matrix[depeak_subtype$Peak.id,]
  subtype_annotation=data.frame(Sample=rownames(subtype_sample_id),subtype=subtype_sample_id$subtype)
  subtype_annotation=merge(subtype_annotation,TCGA_subtype,all.x=T,by="Sample")
  rownames(clinic_data_for_analysis)=paste("X",clinic_data_for_analysis$t_id,sep="")
  subtype_annotation=cbind(subtype_annotation,clinic_data_for_analysis[subtype_annotation$Sample,])
  subtype_annotation=subtype_annotation[,c("subtype","bailey","Collisson","Moffitt","gender","age","fenhua_class","smoke",
                                           "drink","nerve","xueguan","node","TNM_3")]
  
  subtype_peak_annotation=data.frame(Signature_in_subtype=depeak_subtype$Subtype)
  rownames(subtype_peak_annotation)=depeak_subtype$Peak.id
  print("start to draw")
  ann_colors = list(
    TNM_3 = c(`1`='#66C2A5', `2`='#BD0026'),
    node =  c(`0`='#66C2A5', `1`='#BD0026'),
    xueguan = c(`0`='#66C2A5', `1`='#BD0026'),
    nerve = c(`0`='#66C2A5', `1`='#BD0026'),
    drink = c(`0`='#66C2A5', `1`='#BD0026'),
    smoke = c(`0`='#66C2A5', `1`='#BD0026'),
    fenhua_class = c(`1`='#00A2E8', `2`='#13238B',`3`='#FB8072'),
    age = c(`0`='#66C2A5', `1`='#BD0026'),
    gender = c( `1`='#66C2A5',`2`='#BD0026'),
    Moffitt = c(Classical='#FB8072', BasalLike='#3690C0'),
    Collisson = c(Exocrine_like_PDA='#253494', Classical_PDA='#BD0026',QM_PDA='#006837'),
    bailey = c(ADEX='#FF7F27', Immunogenic='#ADD2E5',Pancreatic_Progenitor='#B5E61D',Squamous='#4A9471'),
    subtype = c(s1='#BB3A27', s2='#006AB1',s3='#E18625'),
    Stage= c(I='#BB3A27', II='#006AB1',III='#E18625',IV="#006837"),
    Signature_in_subtype=c(s1='#BB3A27', s2='#006AB1',s3='#E18625',s4="#006837",s5="#B5E61D")
  )                    
  
  de_peak_exp_plot=de_peak_exp_plot[,rownames(subtype_sample_id[order(subtype_sample_id$subtype),])]
  
  pheatmap(log2(de_peak_exp_plot+1),scale='row',border_color = FALSE,
           color = colorRampPalette(c(rep('#0a268a',5),
                                      'white',
                                      rep('#830502',5)))(50),
           annotation_col = subtype_annotation,annotation_row=subtype_peak_annotation,
           annotation_colors = ann_colors,
           show_rownames = F,cluster_cols = F,cluster_rows = F,fontsize = 9,
           filename =  paste(output_dir,"/subtype_heatmap.png",sep=""),width = 12,height = 18)
  return(depeak_subtype)
}


heatmap_for_subtype_for_signature=function(subtype_sample_id,peak_matrix,
                                           TCGA_subtype,clinic_data_for_analysis,output_dir,draw_peak){
  peak_matrix=peak_matrix[draw_peak,]
  jj=1
  number_cluster=unique(subtype_sample_id$subtype)
  depeak_subtype=c()
  while(jj<=length(number_cluster)){
    group_p=rownames(subtype_sample_id[which(subtype_sample_id$subtype==unique(subtype_sample_id$subtype)[jj]),])
    group_non_p=rownames(subtype_sample_id[which(subtype_sample_id$subtype!=unique(subtype_sample_id$subtype)[jj]),])
    #subtype heatmap
    i=1
    n=nrow(peak_matrix)
    depeak_temp=c()
    while(i<=n){
      pvalue=wilcox.test(as.numeric(peak_matrix[i,group_p]),as.numeric(peak_matrix[i,group_non_p]))
      pvalue=pvalue$p.value
      log2FC=log2(mean(as.numeric(peak_matrix[i,group_p]))/mean(as.numeric(peak_matrix[i,group_non_p])))
      depeak_temp=rbind(depeak_temp,c(rownames(peak_matrix)[i],log2FC,pvalue))
      i=i+1
    }
    depeak_temp=data.frame(depeak_temp)
    depeak_temp$X2=as.numeric(depeak_temp$X2)
    depeak_temp$X3=as.numeric(depeak_temp$X3)
    depeak_temp=data.frame(cbind(depeak_temp,FDR=p.adjust(depeak_temp$X3,method="fdr")))
    depeak_temp=depeak_temp[which(depeak_temp$FDR<0.05 ),] #for all signature
    #depeak_temp=depeak_temp[which(depeak_temp$X3<0.05 & abs(depeak_temp$X2)>0.58),] #for 306 signature
    
    if(nrow(depeak_temp)<1){
      return(print("No peak in a subtype!"))
    }
    depeak_temp=depeak_temp[order(depeak_temp$X2,decreasing = T),]
    depeak_subtype=rbind(depeak_subtype,cbind(depeak_temp,Subtype=unique(subtype_sample_id$subtype)[jj]))
    print(jj)
    jj=jj+1
  }
  depeak_subtype=data.frame(depeak_subtype)
  colnames(depeak_subtype)=c("Peak.id","log2FC","P.value","FDR","Subtype")
  depeak_subtype=depeak_subtype[order(abs(depeak_subtype$log2FC),decreasing = T),]
  depeak_subtype=depeak_subtype[!duplicated(depeak_subtype$Peak.id),]
  depeak_subtype=depeak_subtype[order(depeak_subtype$Subtype,depeak_subtype$log2FC),]
  print("order the data")
  de_peak_exp_plot=peak_matrix[depeak_subtype$Peak.id,]
  subtype_annotation=subtype_sample_id
  
  subtype_peak_annotation=data.frame(Signature_in_subtype=depeak_subtype$Subtype)
  rownames(subtype_peak_annotation)=depeak_subtype$Peak.id
  print("start to draw")
  ann_colors = list(
    TNM_3 = c(`0`='#66C2A5', `1`='#BD0026'),
    node =  c(`0`='#66C2A5', `1`='#BD0026'),
    xueguan = c(`0`='#66C2A5', `1`='#BD0026'),
    nerve = c(`0`='#66C2A5', `1`='#BD0026'),
    drink = c(`0`='#66C2A5', `1`='#BD0026'),
    smoke = c(`0`='#66C2A5', `1`='#BD0026'),
    fenhua_class = c(`0`='#00A2E8', `1`='#13238B',`3`='#FB8072'),
    age = c(`0`='#66C2A5', `1`='#BD0026'),
    gender = c( `1`='#66C2A5',`2`='#BD0026'),
    moffitt = c(classic='#FB8072', basal_like='#3690C0'),
    collisson = c(Exocrine_like_PDA='#253494', Classical_PDA='#BD0026',QM_PDA='#006837'),
    bailey = c(ADEX='#FF7F27', Immunogenic='#ADD2E5',Pancreatic_Progenitor='#B5E61D',Squamous='#4A9471'),
    subtype = c(s2='#BB3A27', s1='#006AB1',s3='#1c7911'),
    Stage= c(I='#BB3A27', II='#006AB1',III='#E18625',IV="#006837"),
    Signature_in_subtype=c(s2='#BB3A27', s1='#006AB1',s3='#1c7911',s4="#006837",s5="#B5E61D")
  )                    
  
  de_peak_exp_plot=de_peak_exp_plot[,rownames(subtype_sample_id[order(subtype_sample_id$subtype),])]
  
  pheatmap(log2(de_peak_exp_plot+1),scale='row',border_color = FALSE,
           color = colorRampPalette(c(rep('#0a268a',5),
                                      'white',
                                      rep('#830502',5)))(50),
           annotation_col = subtype_annotation,annotation_row=subtype_peak_annotation,
           annotation_colors = ann_colors,
           show_rownames = F,cluster_cols = F,cluster_rows = F,fontsize = 9,
           filename =  paste(output_dir,"/subtype_heatmap.png",sep=""),width = 12,height = 15)
  return(depeak_subtype)
}

#WER associated with clinic factor
WER_associated_with_clinic=function(clinic_data_for_analysis,WER,gene_matrix,output_dir){
  clinic_factor=c("gender","age","fenhua_class","smoke","drink","nerve",
                  "xueguan","node","TNM_3")
  p=1
  q=length(WER)
  result=c()
  while(p<=q){
    i=1
    n=length(clinic_factor)
    while(i<=n){
      clinic.temp=cbind(clinic_factor=clinic_data_for_analysis[,clinic_factor[i]],
                        Sample=paste("X",clinic_data_for_analysis$t_id,sep=""))
      clinic.temp=data.frame(clinic.temp)
      clinic_factor.temp=unique(clinic.temp$clinic_factor)
      clinic_factor.temp=clinic_factor.temp[order(clinic_factor.temp,decreasing = T)]
      k=length(clinic_factor.temp)
      if(k>2){
        j=1
        loop=c(1:k,1)
        pvalue.for.print=c()
        WER_clinic_boxplot=c()
        while(j<length(loop)){
          print(clinic_factor.temp[loop[j+1]])
          group_p=gene_matrix[WER[p],intersect(colnames(gene_matrix),
                                               clinic.temp[which(clinic.temp$clinic_factor==clinic_factor.temp[loop[j]]),"Sample"])]
          group_non_p=gene_matrix[WER[p],intersect(colnames(gene_matrix),
                                                   clinic.temp[which(clinic.temp$clinic_factor==clinic_factor.temp[loop[j+1]]),"Sample"])]
          pvalue=wilcox.test(as.numeric(group_p),as.numeric(group_non_p))
          pvalue=round(pvalue$p.value,4)
          pvalue.for.print=paste(pvalue.for.print,"\n",clinic_factor.temp[loop[j]]," vs ",clinic_factor.temp[loop[j+1]]," P value: ",pvalue,sep="")
          WER_clinic_boxplot=rbind(WER_clinic_boxplot,
                                   cbind(value=as.numeric(group_non_p),group=clinic_factor.temp[loop[j+1]]),
                                   cbind(value=as.numeric(group_p),group=clinic_factor.temp[loop[j]]))
          Difftype=if(log2(mean(as.numeric(group_p)/mean(as.numeric(group_non_p))))>0){"Greater"}else{"Less"}
          result=rbind(result,c(WER[p],paste(clinic_factor[i]," ",clinic_factor.temp[loop[j]]," vs ",clinic_factor.temp[loop[j+1]],sep=""),pvalue,Difftype))
          j=j+1
        }
        WER_clinic_boxplot=data.frame(WER_clinic_boxplot)
        WER_clinic_boxplot$value=as.numeric(WER_clinic_boxplot$value)
        
        #ggplot(WER_clinic_boxplot)+geom_boxplot(aes(x=WER_clinic_boxplot$group,y=WER_clinic_boxplot$value,
        #                                            fill=WER_clinic_boxplot$group),outlier.alpha = 0)+
        #  theme_classic()+scale_fill_lancet(name=paste(sep="",WER[p],pvalue.for.print))+xlab(WER[p])+
        #  ylab("Gene expression (FPKM)")
        #ggsave(paste(output_dir,WER[p],"_",clinic_factor[i],".png",sep=""),width = 6,height = 5)
        
      } else {
        group_p=gene_matrix[WER[p],intersect(colnames(gene_matrix),
                                             clinic.temp[which(clinic.temp$clinic_factor==clinic_factor.temp[1]),"Sample"])]
        group_non_p=gene_matrix[WER[p],intersect(colnames(gene_matrix),
                                                 clinic.temp[which(clinic.temp$clinic_factor==clinic_factor.temp[2]),"Sample"])]
        pvalue=wilcox.test(as.numeric(group_p),as.numeric(group_non_p))
        pvalue=round(pvalue$p.value,4)
        WER_clinic_boxplot=rbind(cbind(value=as.numeric(group_p),group=clinic_factor.temp[1]),
                                 cbind(value=as.numeric(group_non_p),group=clinic_factor.temp[2]))
        WER_clinic_boxplot=data.frame(WER_clinic_boxplot)
        WER_clinic_boxplot$value=as.numeric(WER_clinic_boxplot$value)
        ggplot(WER_clinic_boxplot)+geom_boxplot(aes(x=WER_clinic_boxplot$group,y=WER_clinic_boxplot$value,
                                                    fill=WER_clinic_boxplot$group),outlier.alpha = 0)+
          theme_classic()+scale_fill_lancet(name=paste(sep="",WER[p],"P value: ",pvalue))+xlab(WER[p])+
          ylab("Gene expression (FPKM)")
        ggsave(paste(output_dir,WER[p],"_",clinic_factor[i],".png",sep=""),width = 6,height = 5)
        Difftype=if(log2(mean(as.numeric(group_p)/mean(as.numeric(group_non_p))))>0){"Greater"}else{"Less"}
        result=rbind(result,c(WER[p],clinic_factor[i],pvalue,Difftype))
      }
      i=i+1
    }
    p=p+1
  }
  result=data.frame(result)
  colnames(result)=c("WER","Clinic.factor","Pvalue","Difftype")
  return(result)
}


#WER corr with peaks
gene_corr_peak=function(peak_for_analysis,gene_for_analysis,peak_matrix,gene_matrix){
  corr_peak=peak_matrix[peak_for_analysis,]
  corr_gene=gene_matrix[gene_for_analysis,]
  write.table(corr_peak,"temp/corr_peak.txt",quote = F,row.names = T,col.names = T,sep="\t")
  write.table(corr_gene,"temp/corr_gene.txt",quote = F,row.names = T,col.names = T,sep="\t")
  i=1
  n=nrow(corr_gene)
  script=c()
  label=abs(rnorm(1))
  while(i<=n){
    script=rbind(script,paste(sep="","Rscript correlation.R temp/corr_peak.txt temp/corr_gene.txt ",
                              rownames(corr_gene)[i]," ",label," temp/"))
    i=i+1
  }
  write.table(script,"run.sh",quote = F,row.names = F,col.names = F)
  system("cat run.sh| xargs -iCMD -P50 bash -c CMD")
  system(paste("cat temp/*",label,"*.txt > temp/temp",sep=""))
  result=read.table("temp/temp",sep="\t",header=F)
  system("rm temp/*")
  colnames(result)=c("Pvalue","Corr","gene","Peak.id")
  return(result)
}

#peak corr gene
peak_corr_gene=function(peak_for_analysis,gene_for_analysis,peak_matrix,gene_matrix){
  corr_peak=peak_matrix[peak_for_analysis,]
  corr_gene=gene_matrix[gene_for_analysis,]
  write.table(corr_peak,"temp/corr_peak.txt",quote = F,row.names = T,col.names = T,sep="\t")
  write.table(corr_gene,"temp/corr_gene.txt",quote = F,row.names = T,col.names = T,sep="\t")
  i=1
  n=nrow(corr_peak)
  script=c()
  label=abs(rnorm(1))
  while(i<=n){
    script=rbind(script,paste(sep="","Rscript correlation.R temp/corr_gene.txt temp/corr_peak.txt ",
                              rownames(corr_peak)[i]," ",label," temp/"))
    i=i+1
  }
  write.table(script,"run.sh",quote = F,row.names = F,col.names = F)
  system("cat run.sh| xargs -iCMD -P50 bash -c CMD")
  system(paste("cat temp/*",label,"*.txt > temp/temp",sep=""))
  result=read.table("temp/temp",sep="\t",header=F)
  #system("rm temp/*")
  colnames(result)=c("Pvalue","Corr","gene","Peak.id")
  return(result)
}


#correlation of peak with host gene
host_gene_corr_peak=function(peak_for_analysis,gene_for_analysis,peak_matrix,gene_matrix){
  peak_for_analysis=cbind(Peak.id=peak_for_analysis,gene=gene_for_analysis)
  peak_for_analysis=na.omit(peak_for_analysis)
  peak_for_analysis=data.frame(peak_for_analysis)
  peak_for_analysis=peak_for_analysis[!duplicated(peak_for_analysis$Peak.id),]
  peak_for_analysis=data.frame(peak_for_analysis)
  corr_peak=peak_matrix[peak_for_analysis$Peak.id,]
  corr_gene=gene_matrix[unique(peak_for_analysis$gene),]
  i=1
  n=nrow(corr_peak)
  result=c()
  while(i<=n){
    
    corr=cor.test(as.numeric(corr_peak[i,]),as.numeric(corr_gene[peak_for_analysis$gene[i],]),method="spearman")
    result=rbind(result,c(rownames(corr_peak)[i],rownames(corr_gene[peak_for_analysis$gene[i],]),
                          corr$p.value,corr$estimate))
    i=i+1
  }
  print(paste(i," gene had calculated",sep=""))
  result=data.frame(result)
  colnames(result)=c("Peak.id","gene","Pvalue","Corr")
  result$Pvalue=as.numeric(result$Pvalue)
  result$Corr=as.numeric(result$Corr)
  return(result)
}


#correlation of WER with host gene
gene_corr_gene=function(gene_for_analysis_1,gene_for_analysis_2,gene_matrix_1,gene_matrix_2){
  corr_gene_1=gene_matrix_1[gene_for_analysis_1,]
  corr_gene_2=gene_matrix_2[gene_for_analysis_2,]
  write.table(corr_gene_1,"temp/corr_gene_1.txt",quote = F,row.names = T,col.names = T,sep="\t")
  write.table(corr_gene_2,"temp/corr_gene_2.txt",quote = F,row.names = T,col.names = T,sep="\t")
  i=1
  n=nrow(corr_gene_1)
  script=c()
  label=abs(rnorm(1))
  while(i<=n){
    script=rbind(script,paste(sep="","Rscript correlation.R temp/corr_gene_2.txt temp/corr_gene_1.txt ",
                              rownames(corr_gene_1)[i]," ",label," temp/"))
    i=i+1
  }
  write.table(script,"run.sh",quote = F,row.names = F,col.names = F)
  system("cat run.sh| xargs -iCMD -P50 bash -c CMD")
  system(paste("cat temp/*",label,"*.txt > temp/temp",sep=""))
  result=read.table("temp/temp",sep="\t",header=F)
  system("rm temp/*")
  colnames(result)=c("Pvalue","Corr","gene","Peak.id")
  return(result)
}

#pathway for cytoscape network
pathway_network=function(david.pathway){
  i=1
  n=nrow(david.pathway)
  result=c()
  while(i<=n){
    j=i+1
    while(j<=n){
      interaction_count=length(intersect(strsplit2(david.pathway[i,'Genes'],split=", "),
                                         strsplit2(david.pathway[j,'Genes'],split=", ")))
      temp=c(david.pathway[i,'Term'],david.pathway[i,'PValue'],david.pathway[j,'Term'],
             david.pathway[j,'PValue'],interaction_count)
      result=rbind(result,temp)
      j=j+1
    }
    i=i+1
  }
  result=data.frame(result)
  colnames(result)=c("pathway_A","Pvalue_A","pathway_B","Pvalue_B","interaction_count")
  result$interaction_count=as.numeric(result$interaction_count)
  return(result)
  
}

#pathway for cytoscape network for clusterprofiler
pathway_network.for.cp=function(david.pathway){
  i=1
  n=length(unique(david.pathway$Description))
  pathway=david.pathway$Description
  result=c()
  while(i<=n){
    j=i+1
    while(j<=n){
      interaction_count=length(intersect(strsplit2(david.pathway[which(david.pathway$Description==pathway[i]),"geneID"],split="[/]"),
                                         strsplit2(david.pathway[which(david.pathway$Description==pathway[j]),"geneID"],split = "[/]")))
      temp=c(pathway[i],david.pathway[which(david.pathway$Description==pathway[i]),"pvalue"][1],
             pathway[j],
             david.pathway[which(david.pathway$Description==pathway[j]),"pvalue"][1]
             ,interaction_count)
      result=rbind(result,temp)
      j=j+1
    }
    print(i)
    i=i+1
  }
  result=data.frame(result)
  colnames(result)=c("pathway_A","Pvalue_A","pathway_B","Pvalue_B","interaction_count")
  result$interaction_count=as.numeric(result$interaction_count)
  return(result)
  
}

#pathway heatmap
pathway_heatmap=function(data,Pathway_split=":",rankby=c(1,2),split.option=FALSE,
                         output_filename,height,width,star.size=5){
  data=data.frame(data)
  pathway_heatmap=data.frame(Term=data$Term,PValue=data$PValue,type=data$type)
  pathway_heatmap$PValue=as.numeric(pathway_heatmap$PValue)
  type=unique(data$type)
  type=type[order(type,decreasing = T)]
  
  i=2
  n=length(type)
  pathway_heatmap_for_plot=pathway_heatmap[which(pathway_heatmap$type==type[1]),]
  while(i<=n){
    pathway_heatmap_for_plot=merge(pathway_heatmap_for_plot,
                                   pathway_heatmap[which(pathway_heatmap$type==type[i]),],
                                   by="Term",all=T)
    i=i+1
  }
  
  i=1
  col=grep("PValue",colnames(pathway_heatmap_for_plot))
  while(i<=nrow(pathway_heatmap_for_plot)){
    pathway_heatmap_for_plot[i,which(is.na(pathway_heatmap_for_plot[i,]))]=1
    if(length(which(as.numeric(pathway_heatmap_for_plot[i,col])<0.05))<1){
      pathway_heatmap_for_plot=pathway_heatmap_for_plot[-i,]
      next
    }
    i=i+1
  }
  
  if(split.option==TRUE){
    rownames(pathway_heatmap_for_plot)=capitalize(strsplit2(pathway_heatmap_for_plot$Term,split = Pathway_split)[,2])
  }else{
    rownames(pathway_heatmap_for_plot)=capitalize(pathway_heatmap_for_plot$Term)
  }
  pathway_heatmap_for_plot$Term=NULL
  
  pathway_heatmap_for_plot=pathway_heatmap_for_plot[order(pathway_heatmap_for_plot[,(rankby[1]-1)*2+1],
                                                          pathway_heatmap_for_plot[,(rankby[2]-1)*2+1],
                                                          decreasing = F),]
  
  p_matrix=pathway_heatmap_for_plot[,grep("PValue",colnames(pathway_heatmap_for_plot))]
  colnames(p_matrix)=type
  pheatmap(p_matrix,scale="none",display_numbers = matrix(ifelse(p_matrix < 0.05, "*", ""),nrow(p_matrix)),
           color = colorRampPalette(c("#cb2815","white"))(30),
           filename = output_filename,width = width,height = height,number_color = "black",
           fontsize_number = star.size,
           show_rownames = T,show_colnames = T,cluster_cols = F,cluster_rows = F,fontsize = 15)
  return(rownames(p_matrix))
}

#pathway fisher
pathway_fisher=function(gene.set,gene.list){
  i=1
  n=nrow(gene.set)
  result=c()
  while(i<=n){
    in.list=(intersect(as.character(unique(gene.set[i,-(1:2)])),
                       unique(gene.list)))
    set.length=length(unique(as.character(gene.set[i,-(1:2)])))
    list.length=length(gene.list)
    pvalue=dhyper(length(in.list),set.length,
                  12404,list.length)
    pvalue.greater=dhyper(length(in.list)+2,set.length,
                          12404,list.length)
    pvalue.less=dhyper(length(in.list)+1,set.length,
                       12404,list.length)
    type=c("greater","less")[which.min(c(pvalue.greater,pvalue.less))]
    
    result=rbind(result,c(gene.set[i,1],pvalue,type,
                          paste(collapse = ",",
                                in.list)))
    i=i+1
  }
  result=data.frame(result)
  colnames(result)=c("Term","PValue","type","Gene")
  return(result)
}

#cluster Profiler
run_clusterProfiler=function(gene_list){
  go.output=c()
  go.output.v2=c()
  kegg.output=c()
  kegg.output.v2=c()
  eg <- bitr(unique(gene_list),
             fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  #go
  go.bp=enrichGO(eg$ENTREZID, 
                 OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 1, 
                 qvalueCutoff = 1)
  
  #kegg
  kegg=enrichKEGG(eg$ENTREZID,organism = "hsa", 
                  pAdjustMethod = 'BH',pvalueCutoff = 1, 
                  qvalueCutoff = 1,keyType = 'kegg')
  
  #barplot(kegg,showCategory=20,drop=T)
  #dotplot(kegg,showCategory=50)
  #input=enrichMap(kegg,n=10)
  #cnetplot(kegg, categorySize="pvalue")
  
  #kegg annotation table
  kegg.output=as.data.frame(kegg)
  if(nrow(kegg.output)>=1){
    kegg.expand=kegg.output[which(kegg.output$Count>1),]
    if(nrow(kegg.expand)>=1){
      i=1
      n=nrow(kegg.expand)
      q=ncol(kegg.expand)
      kegg.expand.temp=c()
      while(i<=n){
        temp=t(matrix(rep(as.character(kegg.expand[i,]),kegg.expand[i,"Count"]),q,kegg.expand[i,"Count"]))
        temp[,8]=strsplit2(temp[1,8],split = "/")
        kegg.expand.temp=rbind(kegg.expand.temp,temp)
        i=i+1
      }
      if(length(kegg.expand.temp)>=1){
        kegg.expand.temp=data.frame(kegg.expand.temp)
        colnames(kegg.expand.temp)=colnames(kegg.output)
        kegg.output=kegg.output[which(kegg.output$Count==1),]
        kegg.output=rbind(kegg.output,kegg.expand.temp)
      }
    }
    colnames(kegg.output)[8]="ENTREZID"
    kegg.output=merge(kegg.output,eg,by="ENTREZID",all.x=T)
    
    i=1
    gene.temp=unique(kegg.output$SYMBOL)
    n=length(gene.temp)
    kegg.output.v2=c()
    while(i<=n){
      temp=kegg.output[which(kegg.output$SYMBOL==gene.temp[i]),]
      temp[1,"Description"]=paste(unique(temp$Description),collapse = ";")
      temp=temp[1,]
      kegg.output.v2=rbind(kegg.output.v2,temp)
      i=i+1
    }
  }
  
  
  
  #go annotation table
  go.output=as.data.frame(go.bp)
  if(nrow(go.output)>=1){
    go.expand=go.output[which(go.output$Count>1),]
    if(nrow(go.expand)>1){
      i=1
      n=nrow(go.expand)
      q=ncol(go.expand)
      go.expand.temp=c()
      while(i<=n){
        temp=t(matrix(rep(as.character(go.expand[i,]),go.expand[i,"Count"]),q,go.expand[i,"Count"]))
        temp[,8]=strsplit2(temp[1,8],split = "/")
        go.expand.temp=rbind(go.expand.temp,temp)
        i=i+1
      }
      if( nrow(go.expand.temp)>1){
        go.expand.temp=data.frame(go.expand.temp)
        colnames(go.expand.temp)=colnames(go.output)
        go.output=go.output[which(go.output$Count==1),]
        go.output=rbind(go.output,go.expand.temp)
      }
    }
    colnames(go.output)[8]="ENTREZID"
    go.output=merge(go.output,eg,by="ENTREZID",all.x=T)
    i=1
    gene.temp=unique(go.output$SYMBOL)
    n=length(gene.temp)
    go.output.v2=c()
    while(i<=n){
      temp=go.output[which(go.output$SYMBOL==gene.temp[i]),]
      temp[1,"Description"]=paste(unique(temp$Description),collapse = ";")
      temp=temp[1,]
      go.output.v2=rbind(go.output.v2,temp)
      i=i+1
    }
  }
  result=list(go.output,go.output.v2,kegg.output,kegg.output.v2)
  return(result)
}


#multi cox scale
multi_cox=function(clinic_data_for_analysis,survival_scatter_plot_value.sig,peak_matrix,
                   offset_cutoff=0.01,specific_hyper_peak,normal_specific_peak){
  survival_scatter_plot_value.sig$HR.lower95=as.numeric(survival_scatter_plot_value.sig$HR.lower95)
  survival_scatter_plot_value.sig$HR.upper95=as.numeric(survival_scatter_plot_value.sig$HR.upper95)
  
  i=1
  n=nrow(survival_scatter_plot_value.sig)
  multi_coxph=c()
  while(i<=n){
    peak.temp=(t(peak_matrix[survival_scatter_plot_value.sig$Peak.id[i],]))
    peak.temp=data.frame(ID=rownames(peak.temp),Peak=as.numeric(peak.temp))
    peak.temp=data.frame(peak.temp,clinic_data_for_analysis[peak.temp$ID,])
    coxph.multi.temp=coxph(Surv(peak.temp$month, peak.temp$status) ~ Peak+gender+age+fenhua_class+smoke+drink+nerve+xueguan+node+TNM_3,
                           data = peak.temp)
    multi_HR=summary(coxph.multi.temp)
    coxph.p.temp=multi_HR$coefficients["Peak","Pr(>|z|)"]
    coxph.HR.temp=multi_HR$coefficients["Peak","exp(coef)"]
    coxph.HR.low95.temp=multi_HR$conf.int["Peak","lower .95"]
    coxph.HR.up95.temp=multi_HR$conf.int["Peak","upper .95"]
    multi_coxph=rbind(multi_coxph,c(survival_scatter_plot_value.sig$Peak.id[i],coxph.p.temp,coxph.HR.temp,
                                    coxph.HR.low95.temp,coxph.HR.up95.temp))
    i=i+1
  }
  multi_coxph=data.frame(multi_coxph)
  colnames(multi_coxph)=c("Peak.id","HR.Pvalue.cox","HR.cox","low95.cox","up95.cox")
  multi_coxph$HR.Pvalue.cox=as.numeric(multi_coxph$HR.Pvalue.cox)
  multi_coxph$HR.cox=as.numeric(multi_coxph$HR.cox)
  multi_coxph$low95.cox=as.numeric(multi_coxph$low95.cox)
  multi_coxph$up95.cox=as.numeric(multi_coxph$up95.cox)
  
  i=1
  n=nrow(multi_coxph)
  multi_coxph.non=c()
  multi_coxph.sig=c()
  while(i<=n){
    if(multi_coxph[i,"HR.Pvalue.cox"]<0.05){
      if(multi_coxph[i,"HR.cox"]>1){
        if(round(multi_coxph[i,"low95.cox"],2)>(1+offset_cutoff)){
          multi_coxph.sig=rbind(multi_coxph.sig,multi_coxph[i,])
        } else {
          multi_coxph.non=rbind(multi_coxph.non,multi_coxph[i,])
        } 
      } else {
        if(round(multi_coxph[i,"up95.cox"],2)<(1-offset_cutoff)){
          multi_coxph.sig=rbind(multi_coxph.sig,multi_coxph[i,])
        } else {
          multi_coxph.non=rbind(multi_coxph.non,multi_coxph[i,])
        }
      }
    } else {
      multi_coxph.non=rbind(multi_coxph.non,multi_coxph[i,])
    }
    i=i+1
  }
  
  table5=rbind(cbind(multi_coxph.sig,type="Sig"),cbind(multi_coxph.non,type="Non"))
  table5=cbind(table5,type_2="Unknown")
  rownames(table5)=table5$Peak.id
  table5[intersect(table5$Peak.id,specific_hyper_peak),"type_2"]="Hyper"
  table5[intersect(table5$Peak.id,normal_specific_peak),"type_2"]="Hypo"
  
  table5=data.frame(cbind(table5$Peak.id,
                          anno[table5$Peak.id,c("Gene","Level.2.gene.type","Gene.site")],
                          survival_scatter_plot_value.sig[table5$Peak.id,c("LogRank.P","HR","HR.Pvalue")],
                          HR.CI=paste(round(survival_scatter_plot_value.sig[table5$Peak.id,"HR.lower95"],2),
                                      "-",
                                      round(survival_scatter_plot_value.sig[table5$Peak.id,"HR.upper95"],2),
                                      sep=""),
                          table5$HR.Pvalue.cox,table5$HR.cox,
                          HR.cox.CI=paste(round(table5$low95.cox,2),"-",round(table5$up95.cox,2),sep=""),
                          type=table5$type,
                          FDR=survival_scatter_plot_value.sig[table5$Peak.id,"FDR"],
                          Bonferroni=survival_scatter_plot_value.sig[table5$Peak.id,"Bonferroni"]
  ))
  colnames(table5)=c("Peak.id","Gene","Gene.type","Gene.site",
                     "LogRank.p","HR","HR.Pvalue","HR.CI",
                     "HR.Pvalue.cox","HR.cox","HR.cox.CI","type","LogRank.FDR","LogRank.Bonferroni")
  table5$LogRank.p=as.numeric(table5$LogRank.p)
  table5$HR=as.numeric(table5$HR)
  table5$HR.Pvalue=as.numeric(table5$HR.Pvalue)
  table5$HR.Pvalue.cox=as.numeric(table5$HR.Pvalue.cox)
  table5$HR.cox=as.numeric(table5$HR.cox)
  table5$LogRank.FDR=as.numeric(table5$LogRank.FDR)
  table5$LogRank.Bonferroni=as.numeric(table5$LogRank.Bonferroni)
  return(table5)
  
}


#make randomforest script
build_RF_script=function(gene_matrix,
                         #clinic_data.for.analysis,
                         target.peak.matrix,target.gene
                         #,Clinic.factor.temp
                         ){
  j=1
  k=nrow(target.peak.matrix)
  #clinic_data.for.analysis=clinic_data.for.analysis[colnames(target.peak.matrix),Clinic.factor.temp]
  clinic_data.for.analysis=data.frame(target.gene=t(gene_matrix[target.gene,colnames(target.peak.matrix)]))
  get_script=c()
  while(j<=k){
    peak.temp=rownames(target.peak.matrix)[j]
    target.temp=target.peak.matrix[j,]
    p=1
    q=nrow(clinic_data.for.analysis)
    data_in_script=c()
    while(p<=q){
      data_in_script=append(data_in_script,paste("[",paste(clinic_data.for.analysis[p,],collapse = ","),"]"))
      p=p+1
    }
    temp.script=paste("sed 's/xingyang_data/[",paste(data_in_script,collapse = ","),
                      "]/g' /data/xingyang/NPC_m6A_zms/py.test/py.script/Construct_model_and_output.py |sed 's/xingyang_target/[",
                      paste(target.temp,collapse = ","),
                      "]/g'|sed 's/feature_list_xingyang/[",paste("\"target.gene.",target.gene,"\"",collapse = ",",sep=""),
                      "]/g'|sed 's/temp_dir_xingyang/",
                      gsub("-","_",gsub(":","_",peak.temp)),"_temp/g' >  /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/",
                      gsub("-","_",gsub(":","_",peak.temp)),"_Construct_model_and_output.py",sep="")
    get_script=append(get_script,temp.script)
    j=j+1
  }
  return(get_script)
}


#Organize RF CDF plot
RF_CDF_data=function(contribution,
                     #Clinic.factor.temp,
                     target.gene){
  contribution_peak=strsplit2(contribution[,1],split = "[{]")[,1]
  contribution_peak=strsplit2(contribution_peak,split = "_|/")
  contribution_peak=contribution_peak[,(grep("Construct",contribution_peak[1,])-1):
                                        (grep("Construct",contribution_peak[1,])+1)]
  contribution_peak=paste(contribution_peak[,1],":",contribution_peak[,2],"-",contribution_peak[,3],sep="")
  
  k=length(target.gene)#length(Clinic.factor.temp)+
  contribution[,1]=strsplit2(contribution[,1],split = "[{]")[,2]
  contribution[,k]=gsub("}","",contribution[,k])
  contribution=as.character(as.matrix(contribution))
  contribution=cbind(contribution,rep(contribution_peak,k))
  contribution=data.frame(contribution)
  
  contribution.target.gene=contribution[grep("target.gene",contribution[,1]),]
  contribution.target.gene=cbind(contribution.target.gene,strsplit2(contribution.target.gene[,1],split=":"))
  contribution.target.gene[,4]=as.numeric(contribution.target.gene[,4])
  contribution.target.gene[,3]=gsub("target.gene.","",contribution.target.gene[,3])
  contribution.target.gene[,3]=gsub(" ","",contribution.target.gene[,3])
  contribution.target.gene=data.frame(contribution.target.gene)
  contribution.target.gene$contribution=NULL
  contribution.target.gene$X2=as.numeric(contribution.target.gene$X2)
  contribution.target.gene[which(contribution.target.gene$X2>1),"X2"]=1
  contribution.target.gene[which(contribution.target.gene$X2<(-1)),"X2"]=-1
  contribution.target.gene=contribution.target.gene[!duplicated(paste(contribution.target.gene[,1],":",contribution.target.gene[,2],sep="")),]
  contribution.target.gene_plot=contribution.target.gene[,2:3]
  return(contribution.target.gene_plot)
}

#Organize RF pvalue plot
RF_P_bar=function(RF_p,Clinic.factor.temp,WER,p_cutoff=0.25){
  k=length(Clinic.factor.temp)+length(WER)
  RF_peak=strsplit2(RF_p[,1],split = "\\[")[,1]
  RF_peak=strsplit2(RF_peak,split = "_|/")
  RF_peak=RF_peak[,(grep("Construct",RF_peak[1,])-1):(grep("Construct",RF_peak[1,])+1)]
  RF_peak=paste(RF_peak[,1],":",RF_peak[,2],"-",RF_peak[,3],sep="")
  RF_p[,1]=strsplit2(RF_p[,1],split = "\\[")[,2]
  RF_p[,k]=gsub("\\]","",RF_p[,k])
  RF_p[,k*2]=gsub("\\]","",RF_p[,k*2])
  
  i=1
  n=grep("WER",RF_p[,1])
  for(i in n){
    temp=RF_p[i,1:k]
    temp=ordered(temp,levels=temp)
    temp=rep(temp,2)
    temp=temp[order(temp)]
    RF_p[i,]=as.character(temp)
  }
  
  RF_p=t(RF_p)
  RF_p_plot=RF_p[grep("\\)",RF_p[,2]),]
  RF_p_plot=rbind(RF_p_plot,RF_peak)
  
  n=seq(2,2*q,by=2)
  q=dim(RF_p_plot)[2]/2
  both_sig=c()
  for(i in n){
    if(length(which(as.numeric(gsub(")","",RF_p_plot[,i]))<p_cutoff))>1){
      both_sig=append(both_sig,RF_p_plot[26,i])
    }
  }
  
  n=seq(1,2*q,by=2)
  q=dim(RF_p_plot)[2]/2
  RF_p_plot_dcast=c()
  for(i in n){
    RF_p_plot_dcast=rbind(RF_p_plot_dcast,cbind(RF_p_plot[1:k,i],RF_p_plot[1:k,i+1],RF_p_plot[k+1,i]))
    
  }
  RF_p_plot_dcast[,2]=gsub("\\)","",RF_p_plot_dcast[,2])
  RF_p_plot_dcast[,1]=gsub(" ","",RF_p_plot_dcast[,1])
  RF_p_plot_dcast[,2]=gsub(" ","",RF_p_plot_dcast[,2])
  RF_p_plot_dcast=data.frame(RF_p_plot_dcast)
  RF_p_plot_dcast$X2=as.numeric(RF_p_plot_dcast$X2)
  RF_p_plot_dcast=data.frame(RF_p_plot_dcast,type="Non")
  RF_p_plot_dcast[which(RF_p_plot_dcast$X2<p_cutoff),"type"]="Sig"
  return(RF_p_plot_dcast)
}

#run randomforest
run_RF=function(
  #clinic_data.for.analysis,Clinic.factor.temp
  target.peak.matrix,gene_matrix,target.gene){
  get_script=build_RF_script(
    #clinic_data.for.analysis=clinic_data.for.analysis,Clinic.factor.temp=Clinic.factor.temp
    target.peak.matrix=target.peak.matrix,
    gene_matrix=gene_matrix,target.gene=target.gene)
  if(length(get_script)<1){return(NULL)}
  write.table(get_script,"/data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/get_script.sh",
              quote = F,row.names = F,col.names = F)
  
  #################################################################
  #sh get_script.sh
  #awk puthon3 chr*.py
  #grep "{" for contribution
  #grep "[" for pvalue
  #################################################################
  system("sh /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/get_script.sh")
  system("ls /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/*_Construct_model_and_output.py|awk \'{print \"python3 \" $1\" > \"$1\".result\"}\'|xargs -iCMD -P60 bash -c CMD")
  system("grep {  /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/*_Construct_model_and_output.py.result > /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/WER.contribution.randomforest.txt")
  #system("grep \"^\\[\" /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/*_Construct_model_and_output.py.result > /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/WER.pvalue.randomforest.txt")
  system("rm -r /data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/*_Construct_model_and_output.py*")
  
  #random forest contribution
  contribution=read.table("/data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/WER.contribution.randomforest.txt",sep=",",fill=T)
  
  contribution.WER_plot=RF_CDF_data(contribution,
                                    #Clinic.factor.temp = Clinic.factor.temp,
                                    target.gene=target.gene)
  
  #random forest pvalue
  # RF_p=read.table("/data/xingyang/NPC_m6A_zms/py.test/py.script/combine_temp/WER.pvalue.randomforest.txt",sep=",",fill=T)
  # 
  # RF_p_plot_dcast=RF_P_bar(RF_p,Clinic.factor.temp,WER,p_cutoff=0.1)
  # 
  #result=list(contribution=contribution.WER_plot,RF_p_plot_dcast=RF_p_plot_dcast)
  result=contribution.WER_plot
  return(result)
}


#site distance
site_distance=function(anno_1,anno_2,gtf.temp){
  i=1
  n=nrow(anno_1)
  result_distribution=c()
  while(i<=n){
    temp.chr=anno_2[which(anno_2[,11]==anno_1[i,11]),]
    if(nrow(temp.chr)<1){
      result_distribution=rbind(result_distribution,
                                c(anno_1[i,4],NA))
      i=i+1
      next
    }
    if(anno_1[i,13]=="intron" | nrow(temp.chr)<2){
      temp.dis=anno_1[i,2]-temp.chr[,2]
      temp.dis=temp.dis[which.min(abs(temp.dis))]
      result_distribution=rbind(result_distribution,
                                c(anno_1[i,4],temp.dis))
      i=i+1
      next
    } else {
      temp.dis=anno_1[i,2]-temp.chr[,2]
      temp.dis=temp.dis[which.min(abs(temp.dis))]
      if(abs(temp.dis)<=2000){
        result_distribution=rbind(result_distribution,
                                  c(anno_1[i,4],temp.dis))
        i=i+1
        next
      } else {
        temp.gtf=gtf.temp[which(gtf.temp$Transcript.id==anno_1[i,8]),]
        temp.gtf=temp.gtf[which(temp.gtf$V3!="transcript" & temp.gtf$V3!="stop_codon" & temp.gtf$V3!="start_codon"),]
        temp.gtf=temp.gtf[order(temp.gtf$V5,decreasing = T),]
        temp.gtf=temp.gtf[order(temp.gtf$V4),]
        temp.gtf=temp.gtf[!duplicated(temp.gtf$V4),]
        temp.gtf=temp.gtf[!duplicated(temp.gtf$V5),]
        if(nrow(temp.gtf)<2){
          result_distribution=rbind(result_distribution,
                                    c(anno_1[i,4],temp.dis))
          i=i+1
          next
        }
        temp.intron=cbind(temp.gtf$V1[-1],temp.gtf$V5[-nrow(temp.gtf)],temp.gtf$V4[-1])
        temp.intron=data.frame(temp.intron)
        temp.intron$X3=as.numeric(temp.intron$X3)
        temp.intron$X2=as.numeric(temp.intron$X2)
        temp.intron$X4=anno_1[i,2]-temp.intron$X2
        temp.intron$X5=anno_1[i,2]-temp.intron$X3
        if(length(which(temp.intron$X4*temp.intron$X5<0))>1){#intron site test
          print(i)
          break
        }
        j=1
        k=nrow(temp.chr)
        temp.CDS=c()
        while(j<=k){
          temp.intron$X6=temp.chr[j,2]-temp.intron$X2
          temp.intron$X7=temp.chr[j,2]-temp.intron$X3
          if(length(which(temp.intron$X6*temp.intron$X7<0))>0){#intron site
            temp.dis=anno_1[i,2]-temp.chr[j,2]
            temp.CDS=rbind(temp.CDS,c(temp.chr[j,4],temp.dis))
          } else { # CDS site
            temp.intron$X8=temp.intron$X3-temp.intron$X2
            temp.intron=temp.intron[order(temp.intron$X6),]
            temp.site=which.min(abs(temp.intron$X5))
            temp.site.2=which.min(abs(temp.intron$X6))
            if(temp.site==temp.site.2){
              temp.dis=temp.chr[j,2]-temp.intron[,2]
              temp.dis=temp.dis[which.min(abs(temp.dis))]
              temp.CDS=rbind(temp.CDS,c(temp.chr[j,4],temp.dis))
              j=j+1
              next
            }
            temp.intron.dis=c(temp.site,temp.site.2)
            temp.intron.dis=temp.intron.dis[order(temp.intron.dis)]
            temp.intron.dis=sum(temp.intron[temp.intron.dis,"X8"])
            temp.dis=anno_1[i,2]-temp.chr[j,2]
            if(temp.dis<0){
              temp.dis=temp.dis+temp.intron.dis
            } else{
              temp.dis=temp.dis-temp.intron.dis
            }
            temp.CDS=rbind(temp.CDS,c(temp.chr[j,4],temp.dis))
          }
          j=j+1
        }
        temp.CDS=data.frame(temp.CDS)
        colnames(temp.CDS)=c("Site.id","Distance")
        temp.CDS$Distance=as.numeric(temp.CDS$Distance)
        result_distribution=rbind(result_distribution,
                                  c(anno_1[i,4],temp.CDS[which.min(abs(temp.CDS$Distance)),"Distance"]))
      }
      
    }
    i=i+1
  }
  result_distribution=as.data.frame(result_distribution)
  result_distribution$V2=as.numeric(result_distribution$V2)
  colnames(result_distribution)=c("Peak.id","Distance")
  return(result_distribution)
}


#Coding annotation with stop codon
coding_anno=function(annotation_table,gtf.temp1){
  colnames(gtf.temp1)[11]="temp.site"
  pre_stop_codon=annotation_table[which(annotation_table$Gene.site=="CDS"|annotation_table$Gene.site=="3UTR"),]
  restore_anno=annotation_table[which(annotation_table$Gene.site!="CDS" & annotation_table$Gene.site!="3UTR"),]
  
  pre_stop_codon=merge(gtf.temp1,pre_stop_codon,by="Transcript.id",all.y=T)
  
  pre_stop_codon$temp.start=pre_stop_codon$V4-pre_stop_codon$Start
  pre_stop_codon$temp.end=pre_stop_codon$V5-pre_stop_codon$Start
  pre_stop_codon[which(pre_stop_codon$temp.start>0),"temp.start"]=1
  pre_stop_codon[which(pre_stop_codon$temp.start<0),"temp.start"]=(-1)
  pre_stop_codon[which(pre_stop_codon$temp.end>0),"temp.end"]=1
  pre_stop_codon[which(pre_stop_codon$temp.end<0),"temp.end"]=(-1)
  pre_stop_codon=pre_stop_codon[which((pre_stop_codon$temp.start*pre_stop_codon$temp.end)<=0),]
  pre_stop_codon$temp.site=pre_stop_codon$Gene.site
  pre_stop_codon[which(pre_stop_codon$Gene.site=="3UTR"),"temp.site"]="UTR"
  pre_stop_codon=pre_stop_codon[which(pre_stop_codon$temp.site==pre_stop_codon$V3),]
  pre_stop_codon$temp.stop=NA
  
  i=1
  n=nrow(pre_stop_codon)
  while(i<=n){
    temp.anno=pre_stop_codon[i,]
    if(temp.anno$V3=="CDS"){
      if(temp.anno$V7=="-"){
        pre_stop_codon[i,"temp.stop"]=temp.anno$Start-temp.anno$V4
      } else {
        pre_stop_codon[i,"temp.stop"]=temp.anno$Start-temp.anno$V5
      }
    } else {
      if(temp.anno$V7=="-"){
        pre_stop_codon[i,"temp.stop"]=temp.anno$Start-temp.anno$V5
      } else {
        pre_stop_codon[i,"temp.stop"]=temp.anno$Start-temp.anno$V4
      }
    }
    i=i+1
  }
  
  restore_anno$temp.stop=NA
  
  stop_codon_anno=rbind(pre_stop_codon[,colnames(restore_anno)],restore_anno)
  stop_codon_anno[which(abs(stop_codon_anno$temp.stop)<100),"Gene.site"]="stop_codon"
  
  stop_codon_regoins_anno=stop_codon_anno[which(stop_codon_anno$Gene.type=="coding"),]
  return(stop_codon_regoins_anno)
}



