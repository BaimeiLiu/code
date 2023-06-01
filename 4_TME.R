range(hnsc.tcga.exp)
###########
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

save(tcga.pcd.ssgsea,file = 'raw_data/TCGA/tcga.pcd.ssgsea.RData')

###########
cancer.immu.genesets=readxl::read_excel('raw_data/cancer-immunity cycle.geneSets.PMID30154154.xlsx')
cancer.immu.genesets=data.frame(cancer.immu.genesets)
head(cancer.immu.genesets)

cancer.immu.genesets.list=split(x=cancer.immu.genesets,f=cancer.immu.genesets$Signature.set)
cancer.immu.genesets.list=sapply(cancer.immu.genesets.list, function(x){subset(x,select='GeneSymbol',drop=TRUE)})

tcga.cancer.immu.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = hnsc.tcga.exp
                                                        ,genelist = cancer.immu.genesets.list)

save(tcga.cancer.immu.ssgsea,file = 'raw_data/immune/tcga.cancer.immu.ssgsea.RData')

################
load('raw_data/TCGA/hnsc.tcga.t.exp.RData')
library('oncoPredict')

GDSC2_Expr = readRDS(file=file.path(dir,'Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"Training Data/GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)

calcPhenotype(trainingExprData = as.matrix(GDSC2_Expr),
              trainingPtype = as.matrix(GDSC2_Res),
              testExprData = as.matrix(hnsc.tcga.t.exp),
              batchCorrect = 'eb',  #   "eb" for ComBat
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData' )
################
tcga.exp=hnsc.tcga.exp
dim(tcga.exp)

library('GSVA')
library(GSEABase)
gmtFile='raw_data/c2.cp.kegg.v7.5.1.symbols.gmt'
c2KEGG <- getGmt(gmtFile,
                 collectionType=BroadCollection(category="c2"),
                 geneIdType=SymbolIdentifier())
tcga.kegg.ssgsea <- gsva(as.matrix(tcga.exp),
                         c2KEGG,
                         method = 'ssgsea',
                         min.sz = 10,
                         max.sz = 500,
                         verbose = TRUE)
save(tcga.kegg.ssgsea,file='raw_data/TCGA/tcga.kegg.ssgsea.RData')

###################
geneSets <- getGmt("./raw_data/h.all.v7.5.1.symbols.gmt")

tcga.exp.h.all <- gsva(as.matrix(tcga.exp), geneSets,min.sz=10,method='ssgsea')
tcga.exp.h.all <- data.frame(tcga.exp.h.all,check.names = F)
save(tcga.exp.h.all,file = 'raw_data/TCGA/tcga.exp.h.all.RData')

##################
library(IOBR)
#### ESTIMATE
tcga.exp.estimate<-deconvo_estimate(eset=tcga.exp)
save(tcga.exp.estimate,file='raw_data/immune/tcga.exp.estimate.RData')
### CIBERSORT
tcga.exp.cibersort<-deconvo_cibersort(eset=tcga.exp,arrays=T)
save(tcga.exp.cibersort,file='raw_data/immune/tcga.exp.cibersort.RData')
############ TIMER 
tcga.exp.timer<-deconvo_timer(eset=as.matrix(tcga.exp),indications=rep('hnsc',ncol(tcga.exp)))
save(tcga.exp.timer,file='raw_data/immune/tcga.exp.timer.RData')
############ EPIC 
tcga.exp.epic<-deconvo_epic(eset=as.matrix(tcga.exp),tumor = TRUE)
save(tcga.exp.epic,file='raw_data/immune/tcga.exp.epic.RData')
############ MCP-counter 
tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(tcga.exp))
save(tcga.exp.mcp,file='raw_data/immune/tcga.exp.mcp.RData')