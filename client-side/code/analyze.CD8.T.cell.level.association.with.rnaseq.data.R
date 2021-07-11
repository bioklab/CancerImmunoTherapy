require(plyr)
require(dplyr)
require(DESeq2)
require(foreach)
require(reshape2)
require(data.table)
require(fgsea)

load('client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')

TCGA.cancer.file.meta <- list(
  ACC  = 'Adrenocortical Cancer.RData',
  BLCA = 'Bladder Urothelial Carcinoma.RData',
  BRCA = 'Breast Invasive Carcinoma.RData',
  CESC = 'Cervical & Endocervical Cancer.RData',
  COAD = 'Colon Adenocarcinoma.RData',
  GBM  = 'Glioblastoma Multiforme.RData',
  HNSC = 'Head & Neck Squamous Cell Carcinoma.RData',
  KICH = 'Kidney Chromophobe.RData',
  KIRC = 'Kidney Clear Cell Carcinoma.RData',
  KIRP = 'Kidney Papillary Cell Carcinoma.RData',
  LGG  = 'Brain Lower Grade Glioma.RData',
  LIHC = 'Liver Hepatocellular Carcinoma.RData',
  LUAD = 'Lung Adenocarcinoma.RData',
  LUSC = 'Lung Squamous Cell Carcinoma.RData',
  OV   = 'Ovarian Serous Cystadenocarcinoma.RData',
  PRAD = 'Prostate Adenocarcinoma.RData',
  READ = 'Rectum Adenocarcinoma.RData',
  SKCM = 'Skin Cutaneous Melanoma.RData',
  THCA = 'Thyroid Carcinoma.RData',
  UCEC = 'Uterine Corpus Endometrioid Carcinoma.RData',
  UCS  = 'Uterine Carcinosarcoma.RData'
)



cancer.type.list <- TCGA.CD8.T.cell.level.df$cancer.type %>% unique %>% as.character

TCGA.rnaseq.analysis.list <- foreach(cancer.type = cancer.type.list) %do% {
    df                <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancery.type == cancer.type,]  
    load(sprintf('server-side//RData//TCGA/%s',TCGA.cancer.file.meta[[cancer.type]]))
  
    data.for.fit   <- TCGA.CD8.T.cell.level.df[TCGA.CD8.T.cell.level.df$cancer.type == cancer.type,]
    median.expr    <- apply(log2.fpkm.matrix,1,median)
    expressed.gene <- names(median.expr)[median.expr > 1]
    expr.matrix    <- log2.fpkm.matrix[expressed.gene,rownames(data.for.fit)]
    data.for.fit$TME.purity <- 1- data.for.fit$tumor.purity
    
    cor.vec <- foreach(g = rownames(expr.matrix),.combine='c') %do% {
        data.for.fit$expr <- expr.matrix[g,]  
        loess.fit         <- loess(data = data.for.fit,formula=expr~tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
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
names(TCGA.rnaseq.analysis.list) <- cancer.type.list

save(file='client-side/output/analyze.CD8.T.cell.level.association.with.rnaseq.data.output/analyze.CD8.T.cell.level.association.with.rnaseq.data.RData',list=c('TCGA.rnaseq.analysis.list'))


df.back <- foreach(item = TCGA.rnaseq.analysis.list,.combine='rbind') %do% {
    item$gene.id <- rownames(item)  
    item
  
}
df <- df.back[df.back$cor.value >=1.5,]
df <- df.back[df.back$cor.value <= -1.5,]

tmp <- ddply(df,.(gene.id), nrow )
tmp <- tmp[order(tmp$V1,decreasing = TRUE),]





pan.cancer.rnaseq.analysis.rs <- tmp[tmp$n >=3,]


x <- TCGA.rnaseq.analysis.list$LGG
x <- TCGA.rnaseq.analysis.list$LIHC

ref: "Tumor-derived CXCL5 promotes human colorectal cancer metastasis through activation of the ERK/Elk-1/Snail and AKT/GSK3β/β-catenin pathways"





