#######
tcga.pcd.cox=cox_batch(t(scale(t(as.matrix(tcga.t.pcd.exp))))
                       ,time = hnsc.tcga.t.exp.os$OS.time
                       , event = hnsc.tcga.t.exp.os$OS)

table(tcga.pcd.cox$p.value<0.05)
table(tcga.pcd.cox$p.value<0.01)
table(tcga.pcd.cox$p.value<0.001)

tcga.pcd.cox=tcga.pcd.cox[order(tcga.pcd.cox$HR,decreasing = T),]
tcga.pcd.cox.sig=tcga.pcd.cox[which(tcga.pcd.cox$p.value<cox.pval),]
nrow(tcga.pcd.cox.sig)


intersect(rownames(tcga.pcd.cox.sig),rownames(tcga.pcd.deg.sig))

############
tcga.pcd.cox.uni=tcga.pcd.cox.sig
tcga.pcd.cox.uni=signif(tcga.pcd.cox.uni,digits=3)
tcga.pcd.cox.uni$`Hazard Ratio(95%CI)`=paste0(tcga.pcd.cox.uni$HR,"(",tcga.pcd.cox.uni$`Low 95%CI`, "-", tcga.pcd.cox.uni$`High 95%CI`, ")")
dim(tcga.pcd.cox.uni)

my_Forestplot(tcga.pcd.cox.uni[rownames(tcga.pcd.cox.uni),c(2:4,1)]
              ,outFile='PDFs/FigS1.pdf',width = 6,height = 6)

my_Forestplot(tcga.pcd.cox.uni[rownames(tcga.pcd.cox.uni),c(2:4,1)]
              ,outFile='02_ConsensusClusterPlus/FigS1.pdf',width = 6,height = 6)

######################
dim(tcga.pcd.cox.sig)
lapply(pcd.genesets.list,function(x){length(intersect(rownames(tcga.pcd.cox.sig),x))})
names(pcd.genesets.list)[lapply(pcd.genesets.list,function(x){length(intersect(rownames(tcga.pcd.cox.sig),x))})>0]

###############
tcga.t.pcd=tcga.t.pcd.exp[rownames(tcga.pcd.cox.sig),]
dim(tcga.t.pcd)
#############################################
################## 分子分型
#############################################
getwd()
setwd('./02_ConsensusClusterPlus')
library(ConsensusClusterPlus)
clusterAlg_use='kmdist'
distance_use='pearson'
#########
clusterAlg_use
distance_use
#######
dim(tcga.t.pcd)
df_exp=tcga.t.pcd
df_exp=as.matrix(df_exp)
df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
df_exp=as.dist(1-cor(df_exp,method = 'spearman'))
dim(df_exp)

clust_subtype = ConsensusClusterPlus(df_exp
                                     , maxK = 10, reps = 500
                                     , pItem = 0.8, pFeature = 1
                                     , title = "TCGA_subtype"
                                     , clusterAlg = clusterAlg_use
                                     , seed = 12345.6789
                                     # , distance = distance_use
                                     # , innerLinkage='ward.D2'
                                     , plot = "pdf", writeTable = T)

k1=4
tcga.subtype=data.frame(clust_subtype[[k1]]$consensusClass)
colnames(tcga.subtype)=c('Cluster')
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
tcga.subtype$Cluster=gsub("CC","C",tcga.subtype$Cluster)
table(tcga.subtype)

library(survcomp)
fig2b=ggplotKMCox(data.frame(time = hnsc.tcga.t.exp.os$OS.time/365
                             , event = hnsc.tcga.t.exp.os$OS
                             , groups=tcga.subtype$Cluster)
                  , title='Cluster'
                  , pal = subtype.color
                  , labs = c('C1','C2','C3','C4')
                  , add_text = '')
fig2b

fig2c=ggplotKMCox(data.frame(time = hnsc.tcga.t.cli$A8_New_Event_Time/365
                             , event = hnsc.tcga.t.cli$A8_New_Event
                             , groups=tcga.subtype$Cluster)
                  , title='Cluster'
                  , pal = subtype.color
                  , labs = c('C1','C2','C3','C4')
                  , add_text = '')
fig2c

library(ggbiplot)
dim(tcga.t.pcd)
tcga.pca<-prcomp(t(tcga.t.pcd), scale=T)
pca.result<-as.data.frame(tcga.pca$x)
pca.result$Subtype<-tcga.subtype$Cluster

fig2a=ggbiplot(tcga.pca, scale=1, groups = tcga.subtype$Cluster,
               ellipse = T,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = subtype.color) +
  # xlim(-3, 3) + ylim(-3, 3) +
  theme_classic() +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  xlab('PCA1') + ylab('PCA2')
fig2a

####################################
library(ComplexHeatmap)
tcga.t.pcd.ssgsea=tcga.t.pcd.ssgsea[,rownames(tcga.subtype)]

column_ha=columnAnnotation(
  foo = anno_block(
    gp = gpar(fill = subtype.color),
    labels = c('C1', 'C2', 'C3','C4'),
    labels_gp = gpar(col = "black", fontsize = 12),
    height = unit(4, "mm")
  ))

T.Stage.color=mycolor[c(1:4)]
N.Stage.color=mycolor[c(1:4)]
M.Stage.color=mycolor[c(1:4)]
Stage.color=mycolor[c(1:4)]
Grade.color=mycolor[c(1:4)]
Age.color=mycolor[c(1,4)]
Gender.color=mycolor[c(1,4)]
Status.color=mycolor[c(1,4)]

names(T.Stage.color)=c('T1','T2','T3','T4')
names(N.Stage.color)=c('N0','N1','N2','N3')
names(M.Stage.color)=c('M0','M1')
names(Stage.color)=c('I','II','III','IV')
names(Grade.color)=c('G1','G2','G3','G4')
names(Age.color)=c('<=60','>60')
names(Gender.color)=c('MALE','FEMALE')
names(Status.color)=c('Alive','Dead')


tcga.t.cli_use$Cluster=tcga.subtype$Cluster
tcga.t.cli_arrange$Cluster=tcga.subtype$Cluster

colnames(tcga.t.cli_use)
colnames(tcga.t.cli_use)[c(3,4,6,8,9,10,12)]

smp.order=dplyr::arrange(tcga.t.cli_use,(Cluster),Stage,Grade
                         ,`T Stage`,`N Stage`,Gender,Age1)
smp.order=rownames(smp.order)

column_ha=HeatmapAnnotation(df = tcga.t.cli_use[smp.order,c(3,4,6,8,9,10,12)]
                            ,col = list(`T Stage`=T.Stage.color
                                        ,`N Stage`=N.Stage.color
                                        ,Stage=Stage.color
                                        ,Grade=Grade.color
                                        ,Gender=Gender.color
                                        ,Age1=Age.color
                                        ,Cluster=subtype.color
                            )
                            , na_col = "grey"
                            , annotation_height = unit(0.01, "mm")
                            , gap = unit(1, 'mm')
)

fig2d=Heatmap(as.matrix(t(scale(t(tcga.t.pcd.ssgsea[,smp.order]))))
              , name = "TPM"
              , cluster_rows = T
              , row_title_gp = gpar(fill = mycolor)
              , show_row_dend = F
              # , column_split = tcga.subtype[colnames(tcga.t.pcd.ssgsea), 'Cluster']
              , cluster_columns = F
              , cluster_column_slices = F
              , show_column_dend = F
              , show_column_names = F
              , col = circlize::colorRamp2(c(-4, 0, 4), c('navy', 'white', 'red'))
              , top_annotation=column_ha
              , border = TRUE)


fig2e=get_PlotMutiBoxplot(t(tcga.t.pcd.ssgsea)[match(rownames(tcga.subtype),colnames(tcga.t.pcd.ssgsea)),]
                          ,group = tcga.subtype
                          ,group.val = 'Cluster'
                          ,ylab = 'ssGSEA score'
                          ,group_cols = subtype.color
                          ,xangle = 45
                          ,legend.pos = 'top')+labs(color='Cluster')
fig2e

fig2abc=mg_merge_plot(fig2a,fig2b,fig2c,nrow = 1,ncol = 3,labels = LETTERS[1:3])

colnames(tcga.t.cli_use)
all(rownames(tcga.t.cli_use)==rownames(tcga.subtype))

tcga.t.cli_use$Cluster=tcga.subtype[rownames(tcga.t.cli_use),]
tcga.t.cli_arrange$Cluster=tcga.subtype[rownames(tcga.t.cli_arrange),]

colnames(tcga.t.cli_use)
table(tcga.t.cli_use$Stage)

colnames(tcga.t.cli_use)
colnames(tcga.t.cli_use)[c(3,4,6,8,9,10,11)]
colnames(tcga.t.cli_arrange)[c(3,4,6,8,9,10,11)]

tcga.subtype.cli.cmp=list()
for(i in c(3,4,6,8,9,10)){
  group.color=mycolor[1:4]
  if(length(na.omit(unique(tcga.t.cli_use[,i])))<=2){
    group.color=mycolor[c(1,3)]
  }
  p=plotMutiBar_tmp(table(tcga.t.cli_use[,i],tcga.t.cli_use$Cluster)
                    ,fill.color = group.color
                    ,isAuto = F,showValue = F
                    ,legTitle=colnames(tcga.t.cli_use)[i])
  tcga.subtype.cli.cmp=c(tcga.subtype.cli.cmp,list(p$Bar))
}
length(tcga.subtype.cli.cmp)

fig2f=mg_merge_plot(tcga.subtype.cli.cmp,nrow = 2,ncol = 3
                    ,labels='F')
fig2f

table(tcga.t.cli_use$Grade,tcga.t.cli_use$Cluster)
