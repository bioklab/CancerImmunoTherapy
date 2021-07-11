require(plyr)
require(dplyr)
require(DESeq2)
require(foreach)
require(reshape2)
require(data.table)
require(fgsea)

load('client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')
load('client-side/output/organize.TCGA.mutation.data.R.output/organize.TCGA.mutation.data.RData')


cancer.type.list <- TCGA.CD8.T.cell.level.df$cancer.type %>% unique %>% as.character

TCGA.mutation.gsea.list <- foreach(cancer.type = cancer.type.list) %do% {
    df                <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancer.type == cancer.type,]  
    #df$TCGA.cohort.id <- sapply(rownames(df),get.TCGA.cohort.id)
    mutation.data     <- TCGA.mutation.data.list[[cancer.type]]
    common.cohort     <- intersect(rownames(df),rownames(mutation.data))
  
    m1 <- match(common.cohort,rownames(df))
    m2 <- match(common.cohort,rownames(mutation.data))
  
    #perform gsea analysis here
    score.for.gsea                   <- df$adjusted.CD8.T.cell.level[m1]
    names(score.for.gsea)            <- common.cohort
    tmp                              <- mutation.data[m2,]
  
    signature.list <- foreach(g= colnames(tmp)) %do% {
        v <- tmp[,g]  
        rownames(tmp)[v==1] 
     }
    names(signature.list) <- colnames(tmp)
    gsea.rs               <- fgsea(stats= score.for.gsea,pathways=signature.list,minSize=5, maxSize=500, nperm=20000)
    gsea.rs               <- gsea.rs[order(gsea.rs$pval),]
    gsea.rs$cancer.type   <- cancer.type
    gsea.rs
}
names(TCGA.mutation.gsea.list) <- cancer.type.list

#Let us try this idea: remove the genes which are NOT expressed in cancer cells


#ref: The Notch/Hes1 Pathway Sustains NF-κB Activation through CYLD Repression in T Cell Leukemia. CYLD gene, hnsc
#ref:The cylindromatosis (CYLD) gene and head and neck tumorigenesis
