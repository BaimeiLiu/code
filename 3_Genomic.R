tcga.tmb=tmb(HNSC,logScale = F)
tcga.tmb=data.frame(tcga.tmb)
rownames(tcga.tmb)=paste0(tcga.tmb$Tumor_Sample_Barcode,"-01")
head(tcga.tmb)

setdiff(rownames(tcga.tmb),rownames(tcga.subtype))
setdiff(rownames(tcga.subtype),rownames(tcga.tmb))
colnames(tcga.tmb)
head(tcga.tmb)

mg_violin_use(data.frame(tcga.subtype[match(rownames(tcga.tmb),rownames(tcga.subtype)),1]
                         ,tcga.tmb[,'total_perMB_log'])
              ,melt = T
              ,ylab = "Tumor mutation load"
              ,jitter=T
              ,group.col = subtype.color
              ,test_method = 'kruskal.test'
              ,cmp_test_method = 'wilcox.test'
              ,legend.pos = 'tr'
              ,show_compare = T)

########
tcga.MATH=math.score(maf = HNSC)
tcga.MATH=data.frame(tcga.MATH)
head(tcga.MATH)
tcga.MATH=data.frame(tcga.MATH)
rownames(tcga.MATH)=paste0(tcga.MATH$Tumor_Sample_Barcode,"-01")
head(tcga.MATH)

mg_violin_use(data.frame(tcga.subtype[match(rownames(tcga.MATH),rownames(tcga.subtype)),1]
                         ,tcga.MATH[,'MATH'])
              ,melt = T
              ,ylab = "Intra-tumor genetic heterogeneity"
              ,jitter=T
              ,group.col = subtype.color
              ,test_method = 'kruskal.test'
              ,cmp_test_method = 'wilcox.test'
              ,legend.pos = 'tr'
              ,show_compare = T)

####################
tcga.immu.lands.p1=data.frame(array=paste0(tcga.immu.lands.p1$`TCGA Participant Barcode`,"-01"),tcga.immu.lands.p1,stringsAsFactors = F,check.names = F)
dim(tcga.immu.lands.p1)
#############
colnames(tcga.immu.lands.p2)[1]=c('array')
dim(tcga.immu.lands.p2)
dim(tcga.immu.lands.p3)

tcga.immu.lands.p2p3=merge(tcga.immu.lands.p2,tcga.immu.lands.p3,by='array',all.x=TRUE,all.y=TRUE,sort = F)
tcga.immu.lands.p2p3=tcga.immu.lands.p2p3[,c('array','LOH_n_seg','LOH_frac_altered','purity','ploidy')]
dim(tcga.immu.lands.p2p3)

setdiff(tcga.immu.lands.p2$array,tcga.immu.lands.p2p3$array)
setdiff(tcga.immu.lands.p3$array,tcga.immu.lands.p2p3$array)

tcga.immu.lands=merge(tcga.immu.lands.p1,tcga.immu.lands.p2p3,by='array',all.x=TRUE,all.y=TRUE,sort = F)
dim(tcga.immu.lands)

setdiff(tcga.immu.lands.p2p3$array,tcga.immu.lands$array)
setdiff(tcga.immu.lands.p1$array,tcga.immu.lands$array)

#####################
table(tcga.immu.lands$`TCGA Study`)

tcga.immu.lands<-tcga.immu.lands[which(tcga.immu.lands$`TCGA Study`=='HNSC'),]
rownames(tcga.immu.lands)=tcga.immu.lands$array
dim(tcga.immu.lands)

table(is.na(match(row.names(tcga.subtype),row.names(tcga.immu.lands))))

tcga.immu.lands=tcga.immu.lands[intersect(row.names(tcga.subtype),row.names(tcga.immu.lands)),]
table(tcga.immu.lands$`TCGA Subtype`)
table(tcga.immu.lands$`Immune Subtype`)

dim(tcga.immu.lands)
colnames(tcga.immu.lands)

col.selected=c('Nonsilent Mutation Rate','SNV Neoantigens','Intratumor Heterogeneity'
               ,'Aneuploidy Score','Number of Segments','Fraction Altered','Homologous Recombination Defects'
               ,'LOH_n_seg','LOH_frac_altered')
tcga.geneAlt.p.all=list()
for(val in col.selected){
  ylab=val
  if(val=='LOH_n_seg'){
    ylab='Numbers of Segs with LOH'
  }else if(val=='LOH_frac_altered'){
    ylab='Fraction of Segs with LOH'
  }else if(val=='Homologous Recombination Defects'){
    ylab='Homologous Recombination Defects'
  }
  p=mg_violin_use(data.frame(tcga.subtype[rownames(tcga.immu.lands),1]
                             ,tcga.immu.lands[,val])
                  ,melt = T
                  ,ylab = ylab
                  ,jitter=F
                  ,group.col = subtype.color
                  ,test_method = 'kruskal.test'
                  ,cmp_test_method = 'wilcox.test'
                  ,legend.pos = 'tr'
                  ,show_compare = T)
  
  plot_ggviolin(data.frame(tcga.subtype[rownames(tcga.immu.lands),1]
                           ,tcga.immu.lands[,val])
                ,melt = T,group.col = subtype.color,ylab=ylab,leg.title = 'Cluster')
  tcga.geneAlt.p.all=c(tcga.geneAlt.p.all,list(p))
}
length(tcga.geneAlt.p.all)

fig3a=mg_merge_plot(tcga.geneAlt.p.all,nrow = 3,ncol = 3
                    ,common.legend = T,labels = LETTERS[1:9])
fig3a
##########
table(tcga.immu.lands$`TCGA Subtype`)
tcga.immu.lands$`TCGA Subtype`=gsub("HNSC.","",tcga.immu.lands$`TCGA Subtype`)

table(tcga.immu.lands$`TCGA Subtype`,tcga.subtype[rownames(tcga.immu.lands),1])
table(tcga.immu.lands$`Immune Subtype`,tcga.subtype[rownames(tcga.immu.lands),1])

fig3b=plotMutiBar_tmp(table(tcga.immu.lands$`TCGA Subtype`
                            ,tcga.subtype[rownames(tcga.immu.lands),1])
                      ,showValue = F
                      ,fill.color = useMyCol("paired",n = 12)[c(3,7,5,8,6)])
fig3b

fig3c=plotMutiBar_tmp(table(tcga.immu.lands$`Immune Subtype`
                            ,tcga.subtype[rownames(tcga.immu.lands),1])
                      ,fill.color = useMyCol("paired",n = 12)[c(3,7,5,8,6)]
                      ,showValue = T
                      ,isAuto = F)
fig3c=fig3c$Bar
fig3c

dataSubt=TCGAbiolinks::TCGAquery_subtype("HNSC")
tcga.subtype.forMut=tcga.subtype
rownames(tcga.subtype.forMut)=substr(rownames(tcga.subtype.forMut),1,12)
writeMatrix(tcga.subtype.forMut,'files/tcga.subtype.forMut.txt')
writeMatrix(tcga.subtype.forMut,'03_Genomic/tcga.subtype.forMut.txt')
