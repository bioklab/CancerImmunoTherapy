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


TCGA.clinical.data.list <- foreach(cancer.type = names(TCGA.cancer.file.meta)) %do% {
    ################### load clinical data ######################################################
    cmd.str              <- sprintf('ls %s/clinical | grep %s',gdac.firehose.folder,cancer.type)
    clinical.data.folder <- system(cmd.str,intern = TRUE)
    cmd.str              <- sprintf('ls %s/clinical/%s | grep picked',gdac.firehose.folder,clinical.data.folder)
    file.list            <- system(cmd.str,intern = TRUE)
    
    clinical.data.df <- foreach(file = file.list,.combine='rbind') %do% {
        file.path    <- sprintf("%s/clinical/%s/%s",gdac.firehose.folder,clinical.data.folder,file)
        data         <- fread(input = file.path) %>% as.data.frame
        data         <- t(data)
        r            <- data[1,]
        data         <- data[2:nrow(data),]
        colnames(data)           <- r
        data                     <- data[,2:ncol(data)]    
        df                       <- data.frame(tcga.cohort.id = toupper(rownames(data)))
        if('years_to_birth' %in% colnames(data)){
            df$years_to_birth        <- as.integer(data[,'years_to_birth'])
        }else{
          df$years_to_birth          <- NA
        }
        df$vital_status          <- as.integer(data[,'vital_status'])
        df$days_to_death         <- as.integer(data[,'days_to_death'])
        df$days_to_last_followup <- as.integer(data[,'days_to_last_followup'])
        df                       <- cbind(df,as.data.frame(data[,5:ncol(data)]))
        rownames(df)             <- df$tcga.cohort.id
        df$tcga.cohort.id        <- NULL
        df
    }
    clinical.data.df$cancer.type <- cancer.type
    clinical.data.df
}
names(TCGA.clinical.data.list) <- names(TCGA.cancer.file.meta)    
  
save(file = 'client-side/output/organize.TCGA.clinical.data.R.output/organize.TCGA.clinical.data.RData',list=c('TCGA.clinical.data.list'))


