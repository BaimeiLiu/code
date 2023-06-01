###################
icg.genes=data.frame(icg.genes)
rownames(icg.genes)=icg.genes$Symbol
head(icg.genes)
setdiff(icg.genes$Symbol,rownames(hnsc.tcga.t.exp))

tcga.t.exp.icg=hnsc.tcga.t.exp[rownames(hnsc.tcga.t.exp) %in% icg.genes$Symbol,]
dim(tcga.t.exp.icg)
top.anno=columnAnnotation(foo = anno_block(
  gp = gpar(fill = subtype.color),
  labels = c('C1','C2','C3','C4'),
  labels_gp = gpar(col = "black", fontsize = 12),
  height = unit(4, "mm")
)) 

fig5a=Heatmap(as.matrix(t(scale(t(tcga.t.exp.icg[,rownames(tcga.subtype)]))))
              , name = "log2(TPM+1)"
              , row_split = icg.genes[rownames(tcga.t.exp.icg), 3]
              , cluster_rows = T
              , cluster_row_slices = T
              , row_title_gp = gpar(fill = mycolor[c(2,3,1)])
              , show_row_dend = F
              , column_split = tcga.subtype$Cluster
              , cluster_columns = F
              , cluster_column_slices = T
              , show_column_dend = F
              , show_column_names = F
              , show_row_names = T
              , row_names_gp = gpar(fontsize = 10)
              , col = circlize::colorRamp2(c(-4, 0, 4), c('navy', 'white', 'red'))
              , top_annotation = top.anno
              , border = TRUE)

ICGs.selected=c('PDCD1','CTLA4','CD274')
setdiff(ICGs.selected,rownames(hnsc.tcga.t.exp))
get_PlotMutiBoxplot(t(tcga.t.exp.icg[ICGs.selected,]),tcga.subtype
                    , legend.pos = 'top'
                    , group_cols = subtype.color
                    , group.val = 'Cluster',ylab = 'log2(TPM+1)',xangle=0)+labs(color='Cluster')

p.all=list()
for(gene in ICGs.selected){
  p=mg_violin_use(data.frame(tcga.subtype,t(tcga.t.exp.icg)[rownames(tcga.subtype),gene])
                  ,test_method = 'wilcox.test'
                  ,cmp_test_method = 'wilcox.test'
                  ,group.col = subtype.color
                  , ylab = gene, melt = T, legend.pos = "none")
  p.all=c(p.all,list(p))
}

fig5b=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'B')
fig5b

############
tcga.icg.exp=tcga.t.exp.icg[ICGs.selected,rownames(tcga.subtype)]
dim(tcga.icg.exp)

tcga.icg.cmp=get_compare_means(data = t(tcga.icg.exp),group = tcga.subtype$Cluster)
tcga.icg.cmp=data.frame(tcga.icg.cmp)
head(tcga.icg.cmp)

tcga.icg.exp_use=tcga.icg.exp
all(rownames(tcga.icg.exp_use)==tcga.icg.cmp$type)
rownames(tcga.icg.exp_use)=paste(tcga.icg.cmp$p.signif,rownames(tcga.icg.exp_use))

tcga.icg.exp_m=data.frame(tcga.subtype,t(tcga.icg.exp_use),check.names = F)
tcga.icg.exp_m=reshape2::melt(tcga.icg.exp_m,id.vars='Cluster')
colnames(tcga.icg.exp_m)=c('Cluster','Gene','Exprs')
head(tcga.icg.exp_m)
dim(tcga.icg.exp_m)

mg_ridges_plot_use(data = tcga.icg.exp_m,melt = T,xlab = 'log2(TPM+1)',col=alpha(subtype.color,alpha = 0.7))


#########################################################
# writeMatrix(t(scale(t(hnsc.tcga.t.exp),scale = F)),outpath = 'raw_data/TCGA/hnsc.tcga.t.exp.zscore.txt')
tcga.tide<-read.csv('raw_data/TCGA-HNSC_TIDE.csv',row.names = 1,stringsAsFactors = F)
tcga.tide=tcga.tide[rownames(tcga.subtype),]
dim(tcga.tide)

tide.selected=c('MDSC','CAF','TAM.M2','Exclusion','Dysfunction','TIDE')

tcga.subtype.tide.p.all=c()
for(i in tide.selected){
  p=plot_ggviolin(data.frame(tcga.subtype$Cluster,tcga.tide[rownames(tcga.subtype),i])
                  ,group.col = subtype.color
                  , ylab = i, melt = T)+labs(fill='Cluster')
  p=mg_violin_use(data.frame(tcga.subtype$Cluster,tcga.tide[rownames(tcga.subtype),i])
                  ,melt = T
                  ,ylab = i
                  ,jitter=T
                  ,group.col = subtype.color
                  ,test_method = 'kruskal.test'
                  ,cmp_test_method = 'wilcox.test'
                  ,legend.pos = 'tr'
                  ,show_compare = T)
  tcga.subtype.tide.p.all=c(tcga.subtype.tide.p.all,list(p))
}
length(tcga.subtype.tide.p.all)
