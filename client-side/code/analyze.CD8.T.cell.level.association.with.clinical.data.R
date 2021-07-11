require(plyr)
require(dplyr)
require(DESeq2)
require(foreach)
require(reshape2)
require(data.table)
require(fgsea)

load('client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')
load('client-side/output/organize.TCGA.clinical.data.R.output/organize.TCGA.clinical.data.RData')

get.TCGA.cohort.id <- function(x) {
    tmp <- strsplit(x,split = '\\-') %>% unlist  
    paste(tmp[1:3],collapse = '-')
  
}

cancer.type.list <- TCGA.CD8.T.cell.level.df$cancer.type %>% unique %>% as.character

pdf('client-side/hehe.pdf')
foreach(cancer.type = cancer.type.list) %do% {
    df                <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancer.type == cancer.type,]  
    df$TCGA.cohort.id <- sapply(rownames(df),get.TCGA.cohort.id)
    clinical.data     <- TCGA.clinical.data.list[[cancer.type]]
    common.cohort     <- intersect(df$TCGA.cohort.id,rownames(clinical.data))
    
    m1 <- match(common.cohort,df$TCGA.cohort.id)
    m2 <- match(common.cohort,rownames(clinical.data))
    
    c.df <- cbind(df[m1,],clinical.data[m2,])
    
    c.df$cancer.type <- NULL
    c.df$TCGA.cohort.id <- NULL
    #perform analysis here
    for(i in 1:ncol(c.df)){
        if(class(c.df[,i]) == 'numeric' | class(c.df[,i]) == 'integer'){    
            plot(x=c.df[,i] %>% factor %>% as.numeric,y=c.df[,'adjusted.CD8.T.cell.level'],xlab=colnames(c.df)[i],main=cancer.type)
        }
        if(class(c.df[,i]) == 'character' ){    
            plot(x=c.df[,i] %>% factor %>% as.numeric,y=c.df[,'adjusted.CD8.T.cell.level'],xlab=colnames(c.df)[i],main=cancer.type)
        }
        if(class(c.df[,i]) == 'factor' & sum( is.na(c.df[,i]) ==TRUE  ) != nrow(c.df) ) {    
            plot(x=c.df[,i] %>% factor %>% as.numeric,y=c.df[,'adjusted.CD8.T.cell.level'],xlab=colnames(c.df)[i],main=cancer.type)
        }
    }
    
}
dev.off()

