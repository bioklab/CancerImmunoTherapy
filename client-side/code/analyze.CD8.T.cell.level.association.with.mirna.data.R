require(plyr)
require(dplyr)
require(DESeq2)
require(foreach)
require(reshape2)
require(data.table)
require(fgsea)

load('client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')
load('client-side/output/organize.TCGA.mirna.data.R.output/organize.TCGA.mirna.data.RData')



cancer.type.list <- TCGA.CD8.T.cell.level.df$cancer.type %>% unique %>% as.character
cancer.type.list <- cancer.type.list[cancer.type.list != 'GBM']

TCGA.mirna.analysis.list <- foreach(cancer.type = cancer.type.list) %do% {
    data.for.fit   <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancer.type == cancer.type,]
    data.for.fit$TME.purity <- 1- data.for.fit$tumor.purity
    
    log2.rpkm.matrix <- TCGA.mirna.data.list[[cancer.type]]
    median.expr      <- apply(log2.rpkm.matrix,1,median)
    expressed.gene   <- names(median.expr)[median.expr > 1]
  
    common.cohort    <- intersect(colnames(log2.rpkm.matrix),rownames(data.for.fit))
    expr.matrix      <- log2.rpkm.matrix[expressed.gene,common.cohort]
    data.for.fit     <- data.for.fit[common.cohort,]
  
    cor.vec <- foreach(g = rownames(expr.matrix),.combine='c') %do% {
        data.for.fit$expr <- 2^expr.matrix[g,]-1 
        loess.fit         <- loess(data = data.for.fit,formula=expr~ tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
        adjusted.expr     <- loess.fit$residuals
        cor(adjusted.expr,data.for.fit$adjusted.CD8.T.cell.level,method='spearman')
    
    #lm.df             <- data.frame(x=adjusted.expr %>% rank,y=data.for.fit$adjusted.CD8.T.cell.score %>% rank)
    #lm.fit.rs         <- lm(data = lm.df,formula =  y ~ x)
    #summary(lm.fit.rs)$coefficients['x',4]
  }
  names(cor.vec) <- rownames(expr.matrix)
  sd.value <- 1.4826 * mad(cor.vec)
  m.value   <- median(cor.vec)
  data.frame(cor.value = (cor.vec - m.value) / sd.value,cancer.type = cancer.type)
  
}
names(TCGA.mirna.analysis.list) <- cancer.type.list

save(file='client-side/output/analyze.CD8.T.cell.level.association.with.mirna.data.output/analyze.CD8.T.cell.level.association.with.mirna.data.RData',list=c('TCGA.mirna.analysis.list'))


df.back <- foreach(item = TCGA.mirna.analysis.list,.combine='rbind') %do% {
  item$gene.id <- rownames(item)  
  item
  
}
df <- df.back[abs(df.back$cor.value) >=1.5,]

df <- df.back[df.back$cor.value <= -1.5,]

tmp <- ddply(df,.(gene.id), nrow )
tmp <- tmp[order(tmp$V1,decreasing = TRUE),]



df <- df.back[abs(df.back$cor.value) >=3,]
"ref: Antitumor immunity is defective in T cell–specific microRNA-155–deficient mice and is rescued by immune checkpoint blockade
"
"ref: miR-142-5p regulates tumor cell PD-L1 expression and enhances anti-tumor immunity.
"

g <- 'hsa-mir-150'
g <- 'hsa-mir-142'
g <- 'hsa-mir-155'
g <- 'hsa-mir-146a'
g  <- 'hsa-mir-342'

tt <- foreach(cancer.type = cancer.type.list,.combine='c') %do% {
  data.for.fit   <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancery.type == cancer.type,]
  data.for.fit$TME.purity <- 1- data.for.fit$tumor.purity
  
  log2.rpkm.matrix <- TCGA.mirna.data.list[[cancer.type]]
  median.expr      <- apply(log2.rpkm.matrix,1,median)
  expressed.gene   <- names(median.expr)[median.expr > 1]
  
  common.cohort    <- intersect(colnames(log2.rpkm.matrix),rownames(data.for.fit))
  expr.matrix      <- log2.rpkm.matrix[expressed.gene,common.cohort]
  data.for.fit     <- data.for.fit[common.cohort,]
  data.for.fit$expr <- expr.matrix[g,]
  median(data.for.fit$expr[data.for.fit$tumor.purity >= 0.9])
}
View(tt)