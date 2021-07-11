require(genefu)
require(plyr)
require(dplyr)
require(GSVA)
require(GSA)
require(cgdsr)
require(ggplot2)
require(foreach)
require(parallel)
require(doParallel)
require(estimate)
require(limma)
require(pheatmap)
source('code/util.R')

load('RData/prepare.meta.data.RData')
load('RData/lincs_signatures_cmpd.RData')



TCGA.RData.file <- system('ls RData| grep RData',wait = TRUE,intern = TRUE)
TCGA.RData.file <- TCGA.RData.file[ (TCGA.RData.file %in% c('signature_meta.RData','prepare.meta.data.RData','lincs_signatures_cmpd.RData')) == FALSE]

DE.rs <- foreach(file = TCGA.RData.file) %do% {
    str <- sprintf('RData/%s',file)
    load(str)
    gender.vec         <- sample.meta.df$gender[match(x=colnames(log2.fpkm.matrix),table = sample.meta.df$sample.id)] %>% factor %>% as.integer
    gender.vec         <- gender.vec -1
    co.variance.matrix <- run.estimate(log2.fpkm.matrix)
    co.variance.matrix <- cbind( co.variance.matrix,gender=gender.vec)
    IFN.gamma.score    <- explore.IFN.gamma.response(log2.fpkm.matrix = log2.fpkm.matrix,co.variance.matrix = co.variance.matrix,geneset.id ='HALLMARK_INTERFERON_GAMMA_RESPONSE')
    
    
    design.df               <- data.frame(StormalScore = co.variance.matrix[,1],ImmuneScore=co.variance.matrix[,2],gender=gender.vec)
    entrez.gene.id          <- ensemble2entrez[rownames(log2.fpkm.matrix)]
    flag                    <- entrez.gene.id %in% rownames(lincs_signatures_cmpd)
    L1000.gene.expr.matrix  <- apply(log2.fpkm.matrix[flag,],2,rank) 
    L1000.gene.expr.matrix  <- L1000.gene.expr.matrix / nrow(L1000.gene.expr.matrix)
    rs.matrix  <- foreach(i=1:nrow(L1000.gene.expr.matrix),.combine='rbind') %do% {
        bias   <- ifelse(min(IFN.gamma.score) > 0,0,1-min(IFN.gamma.score))
        df     <- cbind(design.df,e=L1000.gene.expr.matrix[i,] ,y=(IFN.gamma.score + bias) %>% log10)  
        l1     <- lm(data=df,formula = y ~ StormalScore + ImmuneScore + gender)
        l2     <- lm(data=df,formula = e ~ StormalScore + ImmuneScore + gender)
        lm.fit <- lm(data=df,formula =  y ~ .)
        c(summary(lm.fit)$coefficients['e',c(1,4)],cor(l1$residuals,l2$residuals,method='spearman'))
    }
    rownames(rs.matrix) <- rownames(L1000.gene.expr.matrix)
    colnames(rs.matrix) <- c('effect.size','p.value','adjusted.cor')
    rs.matrix <- cbind(rs.matrix,fdr = p.adjust(rs.matrix[,'p.value'],method='fdr'))
    list(rs.matrix=rs.matrix,sample.no=ncol(log2.fpkm.matrix))
}
names(DE.rs) <- gsub(x=TCGA.RData.file,pattern = '\\.RData',replacement = '')
save(file='RData/pan.cancer.DE.RData',list = c('DE.rs'))
