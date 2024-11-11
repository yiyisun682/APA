library(Seurat)
library(ggplot2)
library(dplyr)
library(ggscatter)

#Load gene expression matrix
load('./APA/RNA_UMI.Rdata')
seurat@meta.data$lineage = as.character(seurat@meta.data$cluster)
Mye = c("Myeloid progenitor","Immature GMP","Early GMP","Intermediate GMP","Late GMP","MDP",'EBM')
Ery = c("MEP","Erythroid progenitor")
Lym = c("Lymphoid progenitor", "MLP","Early BNK","Late BNK")
seurat@meta.data$lineage[seurat@meta.data$lineage %in% Mye] = 'Mye'
seurat@meta.data$lineage[seurat@meta.data$lineage %in% Ery] = 'Ery'
seurat@meta.data$lineage[seurat@meta.data$lineage %in% Lym] = 'Lym'

#Normalize the gene expression matrix
data_exp = data.frame(seurat[["RNA"]]@counts)
norm_data_exp = apply(data_exp, 2, function(x) (x/sum(x))*100000)
tpm = t(as.matrix(log2(norm_data_exp + 1)))
rownames(tpm) = seurat[[]]$lineage
tpm_mean = aggregate(tpm, by = list(row.names(tpm)), FUN = mean, na.rm = TRUE)
rownames(tpm_mean) = tpm_mean$Group.1
tpm_mean = t(tpm_mean[,-1])
colnames(tpm_mean)[2] = 'HSC_MPP'

#Correlation between gene expression and 3'UTR lengths
apa = read.table("./APA/APA_events.txt")
compare = c("Mye_vs_HSC_MPP","Lym_vs_HSC_MPP", "Ery_vs_HSC_MPP")
all = merge(tpm_mean,acc_mean1,by = 'row.names',suffixes = c('_exp','_acc'),all.y=TRUE)
all = na.omit(all)

path = './APA'
for (i in c("Mye", "Lym", "Ery")){
    short = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'shorten'),'ID']
    long = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'lengthen'),'ID']
    all$category = ifelse(all$Row.names %in% long, "lengthen", ifelse(all$Row.names %in% short, "shorten", "non-APA"))
    all$category = factor(all$category, levels = c("lengthen", "shorten", "non-APA"))
    sp = ggscatter(all, x = paste0(i,'_acc'), y = paste0(i,'_exp'), size = 0.8, color = 'category',
        add = "reg.line",  # Add regressin line
        add.params = list(color = "red", fill = "darkgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        palette = c("firebrick3","navy","darkgrey")
    )
    n = sp + stat_cor(method = "pearson", label.x = 0, label.y = 10 ,r.accuracy = 0.01)+
        labs(x='Accesibility',y='Expression',title = paste0(i,' vs HSC/MPP')) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(n,file = paste0(path,i,'_vs_HSC_MPP_exp_ACC_cor.pdf'),width=4, height=4.2)
}

#Differentially expressed genes in differentiated lineages compared with HSC/MPPs
myout = function(x,y){
    marker = FindMarkers(seurat, ident.1 = x, ident.2 = y, min.pct = 0, logfc.threshold = 0)
    marker$compare = paste0(x,'_vs_',y)
    marker$gene = rownames(marker)
    write.csv(marker,file = paste0(path, 'lineage_',x,'_vs_',y,'_DEGs.csv'))
    return(marker)
}
myout("Mye","HSC_MPP")
myout("Lym","HSC_MPP")
myout("Ery","HSC_MPP")

#cdf curves of gene expression level in differentiated lineages
#cdf curves for Mye and Ery
df_p = data.frame()
for (i in c("Mye", "Ery")){
    gene = read.csv(paste0(path,'lineage_', i, '_vs_HSC_MPP_DEGs.csv'), row.names = 1) 
    short = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'shorten'),'ens']
    long = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'lengthen'),'ens']
    all$diff = all[,i] - all$HSC.MPP
    all_merge = merge(gene,all,by.x = 'gene',by.y = 'ID',all.y = TRUE)
    all_merge$category = ifelse(all_merge$ens %in% long, "lengthen", ifelse(all_merge$ens %in% short, "shorten", "non-APA"))
    all_merge$category = factor(all_merge$category, levels = c("lengthen", "shorten", "non-APA"))
    sub = all_merge[which(all_merge$category %in% c("lengthen", "others")),]
    sub = na.omit(sub)
    p_value = ks.test(all_merge[which(all_merge$category == "lengthen"),'avg_log2FC'], all_merge[which(all_merge$category == "others"),'avg_log2FC'])$p.value
    df = data.frame(lineage = i, p_value = p_value)
    df_p = rbind(df, df_p)

    p1 = ggplot(sub, aes(x = avg_log2FC, colour = category)) +  
        stat_ecdf(pad = FALSE)+
        theme_bw() + 
        scale_color_manual(values = c("firebrick3","darkgrey"))+
        theme(axis.text.x = element_text(color='black',hjust = 1,vjust = 1),
	        axis.text.y = element_text(color='black')) +
	    labs(x="avg_log2FC", y="cumulative", title= paste0(i,' vs HSC/MPP CDF')) 
    ggsave(p1,file = paste0(path,'cdf_',i,' vs HSC_MPP_PDUI_exp.pdf'), width=5.5, height=4.2)
}

#cdf curves for Lym
i='Lym'
gene = read.csv(paste0(path,'lineage_',i,'_vs_HSC_MPP_DEGs.csv'),row.names=1)
short = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'shorten'),'ens']
long = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'lengthen'),'ens']
all$diff = all[,i] - all$HSC.MPP
all_merge = merge(gene,all,by.x = 'gene',by.y = 'ID',all.y=TRUE)
all_merge$category <- ifelse(all_merge$ens %in% long, "lengthen", ifelse(all_merge$ens %in% short, "shorten", "non-APA"))
all_merge$category <- factor(all_merge$category, levels = c("lengthen", "shorten", "non-APA"))
print(table(all_merge$category))
sub = all_merge[which(all_merge$category %in% c("shorten", "non-APA")),]
sub= na.omit(sub)
df = data.frame(lineage = i, p_value = p_value)
df_p = rbind(df, df_p)
write.csv(df_p,file = paste0(path,'cdf_',i,' vs HSC_MPP_PDUI_exp_pvalue.csv'))

p1 <- ggplot(sub, aes(x=avg_log2FC, colour = category)) +  
    stat_ecdf(pad = FALSE)+
    theme_bw() + 
    scale_color_manual(values=c("navyblue","darkgrey"))+
    theme(axis.text.x = element_text(color='black',hjust = 1,vjust = 1),
        axis.text.y = element_text(color='black')) +
    labs(x="avg_log2FC",y="cumulative",title= paste0(i,' vs HSC/MPP CDF')) 
ggsave(p1,file = paste0(path,'cdf_',i,' vs HSC_MPP_PDUI_exp.pdf'), width=5.5, height=4.2)
