tcga.cox=cox_batch(t(scale(t(hnsc.tcga.t.exp[tcga.deg.sig,])))
                   ,time = hnsc.tcga.t.exp.os$OS.time
                   ,event = hnsc.tcga.t.exp.os$OS)
dim(tcga.cox)

table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

############
p.cutoff=0.05
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'
table(tcga.cox_use$type)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

##########################
p1 <- ggplot(data = tcga.cox_use,
             aes(x = coef,
                 y = -log10(p.value)))
p1=p1+geom_point(alpha=0.4, size=3.5, aes(color=type))
p1=p1+scale_color_manual(values=c(mg_colors[2],'grey',mg_colors[1]),limits = c("Protective",'None', "Risk"),name='State')
p1=p1+geom_hline(yintercept = -log10(p.cutoff),lty=4,col="black",lwd=0.8)
p1=p1+ylab('-log10(pvalue)')+xlab('Cox coefficient')
# p1=p1+ggrepel::geom_text_repel(data=module.genes.cox_use[which(module.genes.cox_use$p.value<0.05),],aes(label=Gene))
p1=p1+theme_bw()
p1=p1+theme(
  axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
  axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
  legend.text=element_text(face="plain", family="Times", colour="black"  #设置图例的子标题的字体属性
  ),
  legend.title=element_text(face="plain", family="Times", colour="black" #设置图例的总标题的字体属性
  ),
  legend.justification=c(1,1), legend.position=c(1,1)
  ,legend.background = element_rect(fill = NA, colour = NA)
)
p1

table(tcga.cox_use$type)

########################
tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)

library(ggpubr)
ggdotchart(tcga.cox_forVis
           , x="Gene", y="coef",color = "type",
           palette = "aaas",     #配色方案
           legend = "right",     #图例位置
           sorting = "descending",   #上升排序，区别于desc
           add = "segments",    #增加线段
           rotate = TRUE,       #横向显示
           dot.size = tcga.cox_forVis$p.value/2,        #圆圈大小
           # label = round(tcga.cox_use$coef),   #圆圈内数值
           # font.label = list(color="white",size=10, vjust=0.5),   #圆圈内数值字体设置
           ggtheme = theme_pubr())+font("y.text", size = 10, color = "black")+
  theme(legend.position = 'top')+ylab('Cox coefficient')+xlab('')

tcga.exp.sig=hnsc.tcga.t.exp[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)

#################
dim(tcga.exp.sig)
options(ggrepel.max.hnscerlaps = Inf)
tcga.lasso.res=mg_lasso_cox_use(t((t(tcga.exp.sig)))
                                , time = hnsc.tcga.t.exp.os$OS.time
                                , event = hnsc.tcga.t.exp.os$OS
                                , nfolds = 3
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=hnsc.tcga.t.exp[match(tcga.lasso.res$Genes,row.names(hnsc.tcga.t.exp)),]
dim(tcga.exp.for.cox)

lst.modl=createCoxModel_use((t(tcga.exp.for.cox))
                            ,time = hnsc.tcga.t.exp.os$OS.time
                            , event = hnsc.tcga.t.exp.os$OS
                            , isStep = T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

tcga.risk.score=lst.modl$Score
tcga.risk.score=mosaic::zscore(tcga.risk.score)
range(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)
library(dplyr)
fig7c=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  geom_text(aes(label=round(Coef,digits = 3)),color="black",hjust = "left")+
  
  coord_flip() +
  scale_fill_manual(values=pal_nejm(alpha = 0.9)(8)[c(1,4)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Lasso Cox coefficient") +
  theme_classic()+theme(legend.position = c(0,1))

#################
tcga.exp.forCox<- cbind(time=hnsc.tcga.t.exp.os$OS.time,
                        status=hnsc.tcga.t.exp.os$OS,
                        t(hnsc.tcga.t.exp)[rownames(hnsc.tcga.t.exp.os), lst.modl$Genes])
dim(tcga.exp.forCox)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga.exp.forCox))

pdf('PDFs/Fig7C.pdf', height = 4, width = 7,onefile = F)
survminer::ggforest(cox,data=tcga.exp.forCox,noDigits = 3)
dev.off()

############### TCGA
tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=hnsc.tcga.t.exp.os$OS.time/365,
                                                   event=hnsc.tcga.t.exp.os$OS,
                                                   risk=tcga.risk.score), time = "time", event = "event",
                                        variables = c("risk"))
tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
tcga.cutoff=0

fig7D=plotRiskScoreModel_use(riskScore = tcga.risk.score
                             ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                             , time = hnsc.tcga.t.exp.os$OS.time/365
                             , event = hnsc.tcga.t.exp.os$OS
                             , cutoff = tcga.cutoff
                             , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                             , pal = risk.group.color)
fig7d=fig7D[[2]]
fig7d

fig7E=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                             ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                             , time = hnsc.tcga.t.exp.os$OS.time/365
                             , event = hnsc.tcga.t.exp.os$OS
                             , cutoff = tcga.cutoff
                             , labs = c('High','Low')
                             , title = 'RiskType'
                             , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                             , pal = risk.group.color
                             , mks = c(1,3,5))
fig7e=fig7E[[2]]
fig7e

tcga.group=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.group=data.frame(tcga.group)
colnames(tcga.group)='group'
table(tcga.group$group)

############# GSE65858
match(lst.modl$Genes,row.names(gse65858.t.exp))
lst.modl$Genes

row.names(gse65858.t.exp)[grep("ABHD9",row.names(gse65858.t.exp))]
row.names(gse65858.t.exp)[which(row.names(gse65858.t.exp)=="AR")]="AREG"
row.names(gse65858.t.exp)[which(row.names(gse65858.t.exp)=="ABHD9")]="EPHX3"

gse65858.model.dat=gse65858.t.exp[match(lst.modl$Genes,row.names(gse65858.t.exp)),]
gse65858.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                           ,(t(gse65858.model.dat)))
gse65858.risk.score=mosaic::zscore(gse65858.risk.score)

lst.modl$fmla
lst.vd.mod1$fmla

gse65858.cutoff <- survminer::surv_cutpoint(data.frame(time=gse65858.t.cli.os$OS.time/365,
                                                       event=gse65858.t.cli.os$OS,
                                                       risk=gse65858.risk.score), time = "time", event = "event",
                                            variables = c("risk"))
gse65858.cutoff=gse65858.cutoff$cutpoint$cutpoint
gse65858.cutoff=0

fig7F=plotRiskScoreModel_use(riskScore = gse65858.risk.score
                             , dat = t(gse65858.t.exp[intersect(lst.modl$Genes, row.names(gse65858.t.exp)),])
                             , time = gse65858.t.cli.os$OS.time/365
                             , event = gse65858.t.cli.os$OS
                             , cutoff = gse65858.cutoff
                             , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                             , pal = risk.group.color)
fig7f=fig7F[[2]]
fig7f

fig7G=plotCoxModel_Batch_use(riskScore = gse65858.risk.score
                             , dat = t(gse65858.t.exp[intersect(lst.modl$Genes, row.names(gse65858.t.exp)),])
                             , time = gse65858.t.cli.os$OS.time/365
                             , event = gse65858.t.cli.os$OS
                             , cutoff = gse65858.cutoff
                             , labs = c('High','Low')
                             , title = 'RiskType'
                             , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                             , pal = risk.group.color
                             , mks = c(1,3, 5))
fig7g=fig7G[[2]]
fig7g

gse65858.group=ifelse(gse65858.risk.score>gse65858.cutoff,'High','Low')
gse65858.group=data.frame(gse65858.group)
colnames(gse65858.group)='group'
table(gse65858.group)

######################## GSE41613
match(lst.modl$Genes,row.names(gse41613.t.exp))
lst.modl$Genes

rownames(gse41613.t.exp)[grep("C7orf68",rownames(gse41613.t.exp))]
row.names(gse41613.t.exp)[which(row.names(gse41613.t.exp)=="C7orf68")]="HILPDA"

gse41613.model.dat=gse41613.t.exp[intersect(lst.modl$Genes,row.names(gse41613.t.exp)),]

gse41613.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                           ,scale(t(gse41613.model.dat)))
gse41613.risk.score=mosaic::zscore(gse41613.risk.score)
lst.modl$fmla
lst.vd.mod2$fmla

gse41613.cutoff <- survminer::surv_cutpoint(data.frame(time=gse41613.t.cli.os$OS.time/365,
                                                       event=gse41613.t.cli.os$OS,
                                                       risk=gse41613.risk.score), time = "time", event = "event",
                                            variables = c("risk"))
gse41613.cutoff=gse41613.cutoff$cutpoint$cutpoint
gse41613.cutoff=0

fig7H=plotRiskScoreModel_use(riskScore = gse41613.risk.score
                             , dat = t(gse41613.t.exp[intersect(lst.modl$Genes, row.names(gse41613.t.exp)),])
                             , time = gse41613.t.cli.os$OS.time/365
                             , event = gse41613.t.cli.os$OS
                             , cutoff = gse41613.cutoff
                             , pal = risk.group.color)
fig7h=fig7H[[2]]
fig7h

fig7I=plotCoxModel_Batch_use(riskScore = gse41613.risk.score
                             , dat = t(gse41613.t.exp[intersect(lst.modl$Genes, row.names(gse41613.t.exp)),])
                             , time = gse41613.t.cli.os$OS.time/365
                             , event = gse41613.t.cli.os$OS
                             , cutoff = gse41613.cutoff
                             , labs = c('High','Low')
                             , title = 'RiskType'
                             , pal = risk.group.color
                             , mks = c(1, 3, 5))
fig7i=fig7I[[2]]
fig7i

gse41613.group=ifelse(gse41613.risk.score>gse41613.cutoff,'High','Low')
gse41613.group=data.frame(gse41613.group)
colnames(gse41613.group)='group'
table(gse41613.group)

#############
dir.create('08_Model_Clinical')
tcga.t.cli_use$RiskType=tcga.group[rownames(tcga.t.cli_use),'group']
tcga.t.cli_arrange$RiskType=tcga.group[rownames(tcga.t.cli_use),'group']
tcga.t.cli_use$RiskScore=tcga.risk.score
tcga.t.cli_arrange$RiskScore=tcga.risk.score

colnames(tcga.t.cli_arrange)
colnames(tcga.t.cli_arrange)[c(3,4,6,8,9,10,11,12)]

tcga.group.cli.cmp=list()
for(i in c(3:4,6,8,9,10,11,12)){
  group.color=mycolor[1:4]
  if(length(na.omit(unique(tcga.t.cli_arrange[,i])))<=2){
    group.color=mycolor[c(1,4)]
  }
  if(colnames(tcga.t.cli_arrange)[i]=='Cluster'){
    group.color=as.character(subtype.color)
  }
  p1=plotMutiBar_tmp(table(tcga.t.cli_arrange[,i],tcga.t.cli_arrange$RiskType)
                     ,legTitle=colnames(tcga.t.cli_arrange)[i],fill.color = group.color)
  p2=mg_violin_use(data.frame(tcga.t.cli_arrange[,i],tcga.t.cli_arrange$RiskScore)
                   ,test_method = 'kruskal.test'
                   ,cmp_test_method = 'wilcox.test'
                   ,group.col =group.color
                   , xlab=colnames(tcga.t.cli_arrange)[i]
                   , ylab = 'RiskScore', melt = T, legend.pos = NULL)
  tcga.group.cli.cmp=c(tcga.group.cli.cmp,list(p1,p2))
}
length(tcga.group.cli.cmp)

################
range(tcga.risk.score)
RiskScore.color = circlize::colorRamp2(c(-3, 0, 3), c(risk.group.color[2], 'white', risk.group.color[1]))

smp.order=dplyr::arrange(tcga.t.cli_use,RiskScore,RiskType
                         ,desc(Cluster),(Stage),(Grade))
smp.order=rownames(smp.order)

colnames(tcga.t.cli_use)
colnames(tcga.t.cli_use)[c(12,6,8,9,10,13,14)]

column_ha=HeatmapAnnotation(df = tcga.t.cli_use[smp.order,c(12,6,8,9,10,13,14)]
                            ,col = list(Cluster=subtype.color
                                        # ,T.Stage=T.Stage.color
                                        # ,N.Stage=N.Stage.color
                                        ,Stage=Stage.color
                                        ,Grade=Grade.color
                                        ,Gender=Gender.color
                                        ,Age1=Age.color
                                        # ,Status=Status.color
                                        ,RiskType=risk.group.color
                                        ,RiskScore=RiskScore.color
                            )
                            , na_col = "grey"
                            , annotation_height = unit(0.01, "mm")
                            , gap = unit(1, 'mm')
)

fig8b=Heatmap(as.matrix(t(scale(t((hnsc.tcga.t.exp[lst.modl$Genes,smp.order])))))
              , name = "log2(TPM+1)"
              , cluster_rows = T
              , row_title_gp = gpar(fill = mycolor)
              , show_row_dend = F
              , cluster_columns = F
              , show_column_dend = F
              , show_column_names = F
              , col = circlize::colorRamp2(c(-3, 0, 3), c('navy', 'white', 'red'))
              , top_annotation=column_ha
              , border = TRUE)
