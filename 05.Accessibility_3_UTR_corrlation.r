library(ggplot2)
library(dplyr)
library(ggscatter)

#Chromatin accessibility in 3'UTRs of genes
path='./APA/accesibility/'
acc = read.table(paste0(path,'ACC_level.txt'))
meta=read.table(paste0(path,"meta.txt"))
meta$cell = gsub( "R", "", meta$rna_id)
meta_sub = meta[which(meta$cell %in% colnames(acc)), ]
acc = acc[,meta_sub$cell]
info = read.table(paste0(path,'gene_info.txt'))
info$Loci = gsub(":","-",  info$Loci)

#calculating average chromatin accessibility in 3'UTRs of genes
acc = merge(acc,info[,c('ens','Loci')],by.x="row.names", by.y = "Loci", all.x = TRUE)
rownames(acc) = acc$ens
acc = acc[,c(-1,-1100)]
acc = t(acc)
rownames(acc) = meta_sub$lineage
acc_mean = aggregate(acc, by = list(row.names(acc)), FUN = mean, na.rm = TRUE)
rownames(acc_mean)=acc_mean$Group.1
acc_mean = acc_mean[,c(-1)]
acc_mean = data.frame(t(acc_mean))

#PDUI values in 3'UTRs of genes
rownames(info) = info$ens
info = info[,-1]
data = read.table(paste0(path,"scDaPars_imputed_results.txt"))
data = as.matrix(data[,meta$rna_id])
colnames(data) = meta$lineage
dt = apply(data, 1, function(x) tapply(x, colnames(data), mean))
mean = data.frame(t(dt))
mean = merge(mean, info[,c('ens','Loci','ID')],by.x="row.names", by.y = "ens", all.x = TRUE)

apa = read.table(paste0(path,"2lineage_APA_events_filter.txt"))
all = merge(mean,acc_mean,by.x = 'Row.names',by.y = 'row.names',suffixes = c('_pdui','_acc'),all=FALSE)

for (i in c('Lym','Mye','Ery',)){
    apa_sub = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP')),]
    short = apa_sub[which(apa_sub$UTR == 'shorten'),'ens']
    long = apa_sub[which(apa_sub$UTR == 'lengthen'),'ens']
    all$category1 <- ifelse(all$Row.names %in% long, "lengthen", ifelse(all$Row.names %in% short, "shorten", "non-APA"))
    all$category1 = factor(all$category1,levels=c("lengthen","shorten", "non-APA"))

    all$pdui_diff = all[,paste0(i,'_pdui')] - all[,'HSC.MPP_pdui']    
    all$acc_diff = all[,paste0(i,'_acc')] - all[,'HSC.MPP_acc']
    all = all[complete.cases(all[, c("pdui_diff", "acc_diff")]), ]

    sp = ggscatter(all, x = 'acc_diff', y = 'diff_length',
        size = 0.5,
        add = "reg.line",  # Add regressin line
        add.params = list(color = "red", fill = "darkgray"), # Customize reg. line
        fullrange = TRUE, 
        conf.int = FALSE, # Add confidence interval
        color = 'category1',
        palette = c("firebrick3",'navy',"darkgrey")
        )
    n = sp + 
        stat_cor(show.legend = FALSE,method = "pearson", label.x = -0.5, r.accuracy = 0.01, size=5) + 
        labs(x='ΔChromatin accessibility',y = '3UTR length change',title = paste0(i,' vs HSC/MPP')) + 
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(n,file = paste0(path,i,'_ΔACC_Δlength_cor.pdf'),width=4.2,height=4.2)
}

#cdf curves of chromatin accessibility in differentiated lineages
#cdf curves for Mye and Ery
df_p = data.frame()
for (i in c("Mye",  "Ery")){
    all = all[complete.cases(all[,paste0(i,'_acc')], all$HSC.MPP_acc), ]
    all[,'acc_delta'] = all[,paste0(i,'_acc')] - all$HSC.MPP_pdui

    short = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'shorten'),'ens']
    long = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'lengthen'),'ens']
    all$category <- ifelse(all$Row.names %in% long, "lengthen", ifelse(all$Row.names %in% short, "shorten", "non-APA"))
    all$category <- factor(all$category, levels = c("lengthen", "shorten", "non-APA"))

    sub = all[which(all$category %in% c("lengthen", "non-APA")),]
    p_value = ks.test(all[which(all$category == "lengthen"),'acc_delta'], all[which(all$category == "non-APA"),'acc_delta'])$p.value
    df = data.frame(lineage = i, p_value = p_value)
    df_p = rbind(df, df_p)

    p1<- ggplot(sub, aes(x=acc_delta, colour = category)) +  
        stat_ecdf(pad = FALSE)+
        theme_bw() + 
        scale_color_manual(values=c("firebrick3","darkgrey"))+
        theme(axis.text.x = element_text(color='black',hjust = 1,vjust = 1),
	        axis.text.y = element_text(color='black')) +
	    labs(x="acc_delta",y="cumulative",title= paste0(i,' vs HSC/MPP CDF')) 
    ggsave(p1, file = paste0(path,'cdf_delta_',i,' vs HSC_MPP_PDUI_ACC.pdf'),width=5.5,height=4.2)
}

for (i in c("Lym")){
    all <- all[complete.cases(all[,paste0(i,'_acc')], all$HSC.MPP_acc), ]
    all[,'acc_delta'] = all[,paste0(i,'_acc')] - all$HSC.MPP_pdui
    short = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'shorten'),'ens']
    long = apa[which(apa$compare == paste0(i,'_vs_HSC_MPP') & apa$UTR == 'lengthen'),'ens']
    all$category <- ifelse(all$Row.names %in% long, "lengthen", ifelse(all$Row.names %in% short, "shorten", "non-APA"))
    all$category <- factor(all$category, levels = c("lengthen", "shorten", "non-APA"))
    print(table(all$category))
    sub = all[which(all$category %in% c("shorten", "non-APA")),]
    p_value = ks.test(all[which(all$category == "shorten"),'acc_delta'], all[which(all$category == "non-APA"),'acc_delta'])$p.value
    df = data.frame(lineage = i, p_value = p_value)
    df_p = rbind(df, df_p)

    p1<- ggplot(sub, aes(x=acc_delta, colour = category)) +  
        stat_ecdf(pad = FALSE)+
        theme_bw() + 
        scale_color_manual(values=c("navy","darkgrey"))+
        theme(axis.text.x = element_text(color='black',hjust = 1,vjust = 1),
	        axis.text.y = element_text(color='black')) +
	    labs(x="acc_delta",y="cumulative",title= paste0(i,' vs HSC/MPP CDF')) 
    ggsave(p1, file = paste0(path,'cdf_delta_',i,' vs HSC_MPP_PDUI_ACC.png'),width=5.5,height=4.2)
}
