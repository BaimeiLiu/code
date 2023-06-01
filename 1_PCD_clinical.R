pcd.genesets=readxl::read_excel('raw_data/PCD.geneSets.PMID36341760.xlsx')
pcd.genesets=data.frame(pcd.genesets)
pcd.genesets=pcd.genesets[,-1]
head(pcd.genesets)

pcd.genesets.list=list()
pcd.genesets.df=c()
for(i in colnames(pcd.genesets)){
  pcd.genesets.list[[i]]=as.character(na.omit(pcd.genesets[,i]))
  pcd.genesets.df=rbind(pcd.genesets.df,data.frame(PCD=i,Symbol=pcd.genesets[,i],check.names = F,stringsAsFactors = F))
}

tcga.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = hnsc.tcga.exp
                                                ,genelist = pcd.genesets.list)


table(substr(colnames(tcga.pcd.ssgsea),14,15))

######################## PCD在肿瘤样本中的情况
dir.create('01_PCD_Clinical')
##########
table(substr(colnames(tcga.pcd.ssgsea),14,15))
t.smp=colnames(tcga.pcd.ssgsea)[which(substr(colnames(tcga.pcd.ssgsea),14,15)=="01")]
n.smp=colnames(tcga.pcd.ssgsea)[which(substr(colnames(tcga.pcd.ssgsea),14,15)=="11")]
length(t.smp)
length(n.smp)


# p.smp=intersect(substr(t.smp,1,12),substr(n.smp,1,12))
# t.smp=paste0(p.smp,"-01")
# n.smp=paste0(p.smp,"-01")

tcga.smp.group=c(rep("Tumor",length(t.smp)),rep("Normal",length(n.smp)))
tcga.smp.group=data.frame(group=tcga.smp.group)
rownames(tcga.smp.group)=c(t.smp,n.smp)
head(tcga.smp.group)

sampletype.color=alpha(pal_nejm()(8)[c(8,2)],alpha = 0.8)
names(sampletype.color)=c("Tumor","Normal")

library(ComplexHeatmap)
column_ha=HeatmapAnnotation(df = tcga.smp.group
                            ,col = list(group=sampletype.color)
                            , na_col = "grey"
                            , annotation_height = unit(0.01, "mm")
                            , gap = unit(1, 'mm'))
Heatmap(as.matrix(t(scale(t(tcga.pcd.ssgsea)[c(t.smp,n.smp),])))
        , name = "ssGSEA"
        , cluster_rows = T
        , show_row_dend = T
        , column_split = tcga.smp.group$group
        , cluster_columns = T
        , cluster_column_slices = F
        , show_column_dend = F
        , show_column_names = F
        , top_annotation = column_ha
        , col = circlize::colorRamp2(c(-4, 0, 4), c('navy', 'white', 'red'))
        , border = TRUE)

fig1a=mg_PlotMutiBoxplot(t(tcga.pcd.ssgsea)[c(t.smp,n.smp),]
                         , group = tcga.smp.group$group
                         , legend.pos = "top"
                         , add = 'boxplot'
                         , ylab = "ssGSEA score"
                         , group_cols = mycolor[c(1,3)]
                         , test_method = 'wilcox.test')+labs(color = 'Sample Type')
fig1a

####### PCD genes
pcd.genes=unlist(pcd.genesets.list,use.names=FALSE)
length(pcd.genes)
length(unique(pcd.genes))
######
tcga.t.pcd.exp=hnsc.tcga.t.exp[rownames(hnsc.tcga.t.exp) %in% pcd.genes,]
tcga.n.pcd.exp=hnsc.tcga.n.exp[rownames(hnsc.tcga.n.exp) %in% pcd.genes,]
dim(tcga.t.pcd.exp) ## 1215 genes
dim(tcga.n.pcd.exp)

tcga.pcd.deg=mg_limma_DEG_use(exp = cbind(tcga.t.pcd.exp,tcga.n.pcd.exp),group=c(rep("Tumor",ncol(tcga.t.pcd.exp)),rep("Normal",ncol(tcga.n.pcd.exp))),ulab = 'Tumor',dlab = 'Normal')
tcga.pcd.deg$Summary
tcga.pcd.deg.res=tcga.pcd.deg$DEG
tcga.pcd.deg.res=data.frame(tcga.pcd.deg.res,gene=rownames(tcga.pcd.deg.res),Group=rep("Tumor vs Normal",nrow(tcga.pcd.deg.res)))
head(tcga.pcd.deg.res)

fig1b=mg_volcano_custom(diffData = tcga.pcd.deg.res,topGeneN = 10)
fig1b=fig1b+xlab("Sample Type")
fig1b
savePDF('PDFs/Fig1B.pdf',fig1b,height = 4,width = 4)

tcga.pcd.deg.sig=get_DEG(tcga.pcd.deg,logfc.cutoff = log2(2),p.cutoff = 0.05)
dim(tcga.pcd.deg.sig)

table(tcga.pcd.deg.sig$logFC>0)
table(tcga.pcd.deg.sig$logFC<0)

tcga.pcd.deg.sig.enrich=mg_clusterProfiler(genes =rownames(tcga.pcd.deg.sig))

fig1c=barplot(tcga.pcd.deg.sig.enrich$KEGG,width=0.5,showCategory=10)
fig1c=fig1c+theme(axis.text.y=element_text(size=8,face="plain", family="Times", colour="black"))
fig1c

##############################
#################################
tcga.t.pcd.ssgsea=tcga.pcd.ssgsea[,rownames(hnsc.tcga.t.exp.os)]
dim(tcga.t.pcd.ssgsea)

###############
theme_custom_1=theme_classic()+
  theme(axis.text.x= element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(r = 0, l = 0),
        legend.position = c(0.1, 0.2))

theme_custom_2=theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      family = "Times"
    ),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.1, 0.2)
  )

inds=which(!is.na(tcga.t.cli_arrange[,'T.Stage']))
p1=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'T.Stage']
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'T.Stage')

inds=which(!is.na(tcga.t.cli_arrange[,'N.Stage']))
p2=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'N.Stage']
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'N.Stage')
inds=which(!is.na(tcga.t.cli_arrange[,'Stage']))
p3=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'Stage']
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'Stage')
p1=p1+theme_custom_1
p2=p2+theme_custom_1
p3=p3+theme_custom_2

######
fig1d=patchwork::wrap_plots(list(p1,p2,p3)
                            ,plot_annotation(tag_levels = LETTERS[1:3])
                            , ncol = 1, nrow = 3)
fig1d
#############
inds=which(!is.na(tcga.t.cli_arrange[,'Grade']))
p1=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'Grade']
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'Grade')

inds=which(!is.na(tcga.t.cli_arrange[,'Gender']))
p2=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'Gender']
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'Gender')
inds=which(!is.na(tcga.t.cli_arrange[,'Age1']))
p3=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'Age1']
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'Age')

p3=mg_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea[,inds])
                      , group = tcga.t.cli_arrange[inds,'Age1']
                      , legend.pos = ''
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , group_cols = mycolor[1:4]
                      , test_method = 'wilcox.test')+labs(color = 'Age')
p1=p1+theme_custom_1
p2=p2+theme_custom_1
p3=p3+theme_custom_2
######
fig1e=patchwork::wrap_plots(list(p1,p2,p3)
                            ,plot_annotation(tag_levels = LETTERS[1:3])
                            , ncol = 1, nrow = 3)
fig1e

fig1de=mg_merge_plot(fig1d,fig1e,nrow = 1,ncol = 2)
fig1abc=mg_merge_plot(fig1a,fig1b,fig1c,nrow = 3,ncol = 1,labels = LETTERS[1:3])

###########
tcga.t.pcd.ssgsea=tcga.pcd.ssgsea[,rownames(hnsc.tcga.t.exp.os)]
dim(tcga.t.pcd.ssgsea)
tcga.pcd.cox=cox_batch(t(scale(t(as.matrix(tcga.t.pcd.ssgsea))))
                       ,time = hnsc.tcga.t.exp.os$OS.time
                       , event = hnsc.tcga.t.exp.os$OS)

table(tcga.pcd.cox$p.value<0.05)
table(tcga.pcd.cox$p.value<0.01)
table(tcga.pcd.cox$p.value<0.001)

