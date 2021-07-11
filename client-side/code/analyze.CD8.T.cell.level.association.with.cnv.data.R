require(plyr)
require(dplyr)
require(DESeq2)
require(foreach)
require(reshape2)
require(data.table)
require(fgsea)

load('client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')
load('client-side/output/organize.TCGA.cnv.data.R.output/organize.TCGA.cnv.data.RData')


cancer.type.list <- TCGA.CD8.T.cell.level.df$cancer.type %>% unique %>% as.character

get.TCGA.sample.short.id <- function(x) {
  tmp <- strsplit(x,split = '\\-') %>% unlist  
  tmp[4] <- gsub(x=tmp[4],pattern = 'A',replacement = '')
  paste(tmp[1:4],collapse = '-')
  
}


TCGA.cnv.gsea.list <- foreach(cancer.type = cancer.type.list) %do% {
    df                 <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancer.type == cancer.type,]  
    cnv.data           <- TCGA.cnv.data.list[[cancer.type]]
    colnames(cnv.data) <- sapply(colnames(cnv.data),get.TCGA.sample.short.id)
    cnv.data           <- t(cnv.data)
    common.cohort      <- intersect(rownames(df),rownames(cnv.data))
  
    m1 <- match(common.cohort,rownames(df))
    m2 <- match(common.cohort,rownames(cnv.data))
  
    #perform gsea analysis here
    score.for.gsea                   <- df$adjusted.CD8.T.cell.level[m1]
    names(score.for.gsea)            <- common.cohort
    tmp                              <- cnv.data[m2,]
    sd.value                         <- apply(tmp,2,sd)
    tmp                              <- tmp[,sd.value > 0]
    
    rs.df <- foreach(g= colnames(tmp),.combine='rbind') %do% {
        lm.fit.df <- data.frame(y=score.for.gsea,x=tmp[,g])
        #lm.fit.df <- lm.fit.df[abs(lm.fit.df$x) < 2,]
        lm.rs <- lm(data = lm.fit.df,formula = y ~ x)
        #summary(lm.rs)$coefficients[2,c(1,4)] %>% as.data.frame
        data.frame(p.value = summary(lm.rs)$coefficients[2,c(4)], slop = summary(lm.rs)$coefficients[2,c(1)])
    }
    rownames(rs.df) <- colnames(tmp)
    #fdr.vec <- p.adjust(p.value.vec,method='bonferroni')
    rs.df
#     gsea.rs               <- fgsea(stats= score.for.gsea,pathways=signature.list,minSize=5, maxSize=500, nperm=10000)
#     gsea.rs               <- gsea.rs[order(gsea.rs$pval),]
#     gsea.rs$cancer.type   <- cancer.type
#     gsea.rs
}
names(TCGA.cnv.gsea.list) <- cancer.type.list

