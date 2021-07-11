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
  UCS  = 'Uterine Carcinosarcoma.RData'
)

gdac.firehose.folder      <- 'client-side/external.data/from.gdac.firehose'


TCGA.cnv.data.list <- foreach(cancer.type = names(TCGA.cancer.file.meta)) %do% {
  ################### load cnv data ######################################################
    cmd.str              <- sprintf('ls %s/cnv | grep %s',gdac.firehose.folder,cancer.type)
    cnv.data.folder      <- system(cmd.str,intern = TRUE)
    cmd.str              <- sprintf('ls %s/cnv/%s | grep all_lesions.conf_99',gdac.firehose.folder,cnv.data.folder)
    file.list            <- system(cmd.str,intern = TRUE)
    file                 <- file.list[1]
    file.path            <- sprintf("%s/cnv/%s/%s",gdac.firehose.folder,cnv.data.folder,file)
    data                 <- fread(input = file.path) %>% as.data.frame
    flag                 <- grepl(x=colnames(data),pattern='TCGA')
    cnv.data             <- data[,flag]
    flag                 <- grepl(x=data[,'Unique Name'],pattern='CN values')
    cnv.data             <- cnv.data[!flag,]
    rownames(cnv.data)   <- paste(data[!flag,'Unique Name'],data[!flag,'Wide Peak Limits'],sep = "@")
    cnv.data <- as.matrix(cnv.data)
    cnv.data
  
  
  
#     data[,'Locus ID']    <- NULL  
#     data[,'Cytoband']    <- NULL
#     rownames(data)       <- data[,'Gene Symbol']
#     data[,'Gene Symbol'] <- NULL
#     data                 <- as.matrix(data)
#     data
   
}
names(TCGA.cnv.data.list) <- names(TCGA.cancer.file.meta)    

save(file = 'client-side/output/organize.TCGA.cnv.data.R.output/organize.TCGA.cnv.data.RData',list=c('TCGA.cnv.data.list'))


