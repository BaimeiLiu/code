tcga.C1.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype$Cluster
                               ,ulab='C1',dlab = NULL)
tcga.C2.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype$Cluster
                               ,ulab='C2',dlab = NULL)
tcga.C3.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype$Cluster
                               ,ulab='C3',dlab = NULL)
tcga.C4.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype$Cluster
                               ,ulab='C4',dlab = NULL)

kegmt<-read.gmt("origin_datas/h.all.v7.5.1.entrez.gmt") 

tcga.C1.gsea<-GSEA(tcga.C1.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA分析
tcga.C2.gsea<-GSEA(tcga.C2.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA分析
tcga.C3.gsea<-GSEA(tcga.C3.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA分析
tcga.C4.gsea<-GSEA(tcga.C4.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA分析

tcga.C1.gsea.res=tcga.C1.gsea@result
tcga.C2.gsea.res=tcga.C2.gsea@result
tcga.C3.gsea.res=tcga.C3.gsea@result
tcga.C4.gsea.res=tcga.C4.gsea@result
rownames(tcga.C1.gsea.res)


tcga.hallmark.union=Reduce(union,list(C1=rownames(tcga.C1.gsea.res)
                                      , C2 = rownames(tcga.C2.gsea.res)
                                      , C3 = rownames(tcga.C3.gsea.res)
                                      , C4 = rownames(tcga.C4.gsea.res)))


tcga.gsea.heatmap.dat=matrix(c(0),nrow = 4,ncol = length(tcga.hallmark.union))
rownames(tcga.gsea.heatmap.dat)=c('C1', 'C2','C3','C4')
colnames(tcga.gsea.heatmap.dat)=tcga.hallmark.union

tcga.gsea.heatmap.dat[1,match(rownames(tcga.C1.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C1.gsea.res$NES
tcga.gsea.heatmap.dat[2,match(rownames(tcga.C2.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C2.gsea.res$NES
tcga.gsea.heatmap.dat[3,match(rownames(tcga.C3.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C3.gsea.res$NES
tcga.gsea.heatmap.dat[4,match(rownames(tcga.C4.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C4.gsea.res$NES

range(tcga.gsea.heatmap.dat)
fig6a=Heatmap(as.matrix(t(scale(tcga.gsea.heatmap.dat)))
              , name = "NES"
              , rect_gp = gpar(col = "white")
              , cluster_rows = T
              , show_row_dend = F
              , cluster_columns = F
              , show_column_dend = F
              , show_column_names = T
              , show_row_names = T
              , row_names_gp = gpar(fontsize = 8)
              , row_names_side  = c("left")
              , col = circlize::colorRamp2(c(-3, 0, 3), c('navy', 'white', 'red'))
              , border = TRUE)
fig6a
########### 差异表达分析
dim(hnsc.tcga.t.exp)
tcga.t.exp_use=hnsc.tcga.t.exp[rownames(hnsc.tcga.t.exp) %in% ann.pcg$gene_name,]
dim(tcga.t.exp_use)

################
tcga.deg.c1.res=tcga.deg.c1$DEG
tcga.deg.c2.res=tcga.deg.c2$DEG
tcga.deg.c3.res=tcga.deg.c3$DEG
tcga.deg.c4.res=tcga.deg.c4$DEG

tcga.deg.c1.res=data.frame(tcga.deg.c1.res,gene=rownames(tcga.deg.c1.res),Group=rep("C1 vs Other",nrow(tcga.deg.c1.res)))
tcga.deg.c2.res=data.frame(tcga.deg.c2.res,gene=rownames(tcga.deg.c2.res),Group=rep("C2 vs Other",nrow(tcga.deg.c2.res)))
tcga.deg.c3.res=data.frame(tcga.deg.c3.res,gene=rownames(tcga.deg.c3.res),Group=rep("C3 vs Other",nrow(tcga.deg.c3.res)))
tcga.deg.c4.res=data.frame(tcga.deg.c4.res,gene=rownames(tcga.deg.c4.res),Group=rep("C4 vs Other",nrow(tcga.deg.c4.res)))

tcga.deg.res=rbind(tcga.deg.c1.res,tcga.deg.c2.res,tcga.deg.c3.res,tcga.deg.c4.res)

fig6b=mg_volcano_custom(diffData = tcga.deg.res
                        ,tile.col = as.character(subtype.color)
                        ,log2FC.cutoff = log2(1.5),topGeneN = 5)
fig6b

############
tcga.deg.c1$Summary
tcga.deg.c2$Summary
tcga.deg.c3$Summary
tcga.deg.c4$Summary

tcga.deg.c1.sig=get_DEG(tcga.deg.c1,logfc.cutoff=log2(1.5),p.cutoff = 0.05)
tcga.deg.c2.sig=get_DEG(tcga.deg.c2,logfc.cutoff=log2(1.5),p.cutoff = 0.05)
tcga.deg.c3.sig=get_DEG(tcga.deg.c3,logfc.cutoff=log2(1.5),p.cutoff = 0.05)
tcga.deg.c4.sig=get_DEG(tcga.deg.c4,logfc.cutoff=log2(1.5),p.cutoff = 0.05)

#################################
Reduce(intersect,list(rownames(tcga.deg.c1.sig)
                      , rownames(tcga.deg.c2.sig)
                      , rownames(tcga.deg.c3.sig)
                      , rownames(tcga.deg.c4.sig)))
####################
tcga.deg.sig=Reduce(intersect,list(rownames(tcga.deg.c1.sig)
                                   ,rownames(tcga.deg.c2.sig)
                                   ,rownames(tcga.deg.c3.sig)
                                   ,rownames(tcga.deg.c4.sig)))

length(tcga.deg.sig)

writeMatrix(tcga.deg.sig,'files/文件/tcga.deg.sig.txt')
writeMatrix(tcga.deg.sig,'06_Pathway/tcga.deg.sig.txt')
