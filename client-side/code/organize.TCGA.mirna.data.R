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


TCGA.mirna.data.list <- foreach(cancer.type = names(TCGA.cancer.file.meta)) %do% {
  ################### load mirna data ######################################################
  cmd.str              <- sprintf('ls %s/mirna | grep %s',gdac.firehose.folder,cancer.type)
  mirna.data.folder    <- system(cmd.str,intern = TRUE)
  cmd.str              <- sprintf('ls %s/mirna/%s | grep miRseq_RPKM.txt',gdac.firehose.folder,mirna.data.folder)
  file.list            <- system(cmd.str,intern = TRUE)
  file                 <- file.list[1]
  if(FALSE == is.na(file)){
      file.path            <- sprintf("%s/mirna/%s/%s",gdac.firehose.folder,mirna.data.folder,file)
      data                 <- fread(input = file.path) %>% as.data.frame
      rownames(data)           <- data[,'HYBRIDIZATION R']
      data[,'HYBRIDIZATION R'] <- NULL
      data <- as.matrix(data)
      data <- log2(data + 1)
      data
  }else{
      NA  
  }
  
}

names(TCGA.mirna.data.list) <- names(TCGA.cancer.file.meta)  
flag                        <- names(TCGA.mirna.data.list) == 'GBM'
TCGA.mirna.data.list        <- TCGA.mirna.data.list[!flag] # remove GBM, no miRNA seq data from firehose
save(file = 'client-side/output/organize.TCGA.mirna.data.R.output/organize.TCGA.mirna.data.RData',list=c('TCGA.mirna.data.list'))


