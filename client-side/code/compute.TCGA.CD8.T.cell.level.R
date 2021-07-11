require(plyr)
require(dplyr)
require(DESeq2)
require(foreach)
require(reshape2)
require(data.table)
require(fgsea)
#load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')

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
  UCS  = 'Uterine Carcinosarcoma.RData',
  PRAD = 'Pancreatic Adenocarcinoma.RData'
)

gdac.firehose.folder      <- 'client-side/external.data/from.gdac.firehose'
tumor.purity.df           <- read.csv(file = 'client-side/external.data/expired/TCGA.tumor.purity.csv',header=TRUE)
tumor.purity.df$Sample.ID <- gsub(x = tumor.purity.df$Sample.ID,pattern = '01A$',replacement = '01')
rownames(tumor.purity.df) <- tumor.purity.df$Sample.ID

CD8A <- 'ENSG00000153563'
CD8B <- 'ENSG00000172116'
GZMA <- 'ENSG00000145649'
GZMB <- 'ENSG00000100453'
PRF1 <- 'ENSG00000180644'
IFNG <- 'ENSG00000111537'
TRBC2 <- 'ENSG00000211772'
TRBC1 <- 'ENSG00000211751'
NKG2D <- 'ENSG00000213809'

tmp <- foreach(cancer.type = names(TCGA.cancer.file.meta),.combine='rbind') %do% {
    ################### load rna-seq data ##################################################################
    rna.seq.file <- sprintf('server-side//RData/TCGA/%s',TCGA.cancer.file.meta[[cancer.type]])  
    load(rna.seq.file)
  
    ################### compute adjusted CD8 T cell inflitration level######################################
    #CD8.T.cell.level   <- apply(log2.fpkm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)
    CD8.T.cell.level   <- log2.fpkm.matrix[TRBC2,]
    #CD8.T.cell.level   <- log2.fpkm.matrix[IFNG,]
    #CD8.T.cell.level   <- log2.fpkm.matrix[NKG2D,]
    
    
    meta.df              <- data.frame(tumor.purity=tumor.purity.df[names(CD8.T.cell.level),'CPE'],CD8.T.cell.level=CD8.T.cell.level)
    rownames(meta.df)    <- names(CD8.T.cell.level)
    #meta.df$tumor.purity <- log2.fpkm.matrix[TRBC2,rownames(meta.df)] # well, let us try this
    meta.df              <- meta.df[complete.cases(meta.df),]
    loess.fit            <- loess(data = meta.df,formula=CD8.T.cell.level~tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
    sd.value                          <- 1.4826 * mad(loess.fit$residuals)
    md.value                          <- median(loess.fit$residuals)
    meta.df$adjusted.CD8.T.cell.level <- (loess.fit$residuals - md.value)/sd.value
    meta.df$cancer.type              <- cancer.type
    meta.df
  
}

TCGA.CD8.T.cell.level.df <- tmp

save(file = 'client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData',list =c('TCGA.CD8.T.cell.level.df') )
