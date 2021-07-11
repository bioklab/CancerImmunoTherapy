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


TCGA.mutation.data.list <- foreach(cancer.type = names(TCGA.cancer.file.meta)) %do% {
    ################### load mutation data ######################################################
    cmd.str              <- sprintf('ls %s/mutation | grep %s.Mutation',gdac.firehose.folder,cancer.type)
    mutation.data.folder <- system(cmd.str,intern = TRUE)
    cmd.str              <- sprintf('ls %s/mutation/%s | grep maf',gdac.firehose.folder,mutation.data.folder)
    file.list            <- system(cmd.str,intern = TRUE)
    mutation.data.df <- foreach(file = file.list,.combine='rbind') %do% {
        tmp     <- strsplit(file,split = '\\.') %>% unlist 
        tcga.id <- tmp[1]   
        file    <- sprintf("%s/mutation/%s/%s",gdac.firehose.folder,mutation.data.folder,file)
        data    <- fread(input = file) %>% as.data.frame
        data    <- data[data$Variant_Classification != 'Silent',]
        if(nrow(data) == 0){
            return(NULL)
    }
    data$tcga.id <- tcga.id
    data[,c('tcga.id','Hugo_Symbol')]
  }
    mutation.data.df$value           <- 1
    x                                <- dcast(mutation.data.df,formula = tcga.id ~ Hugo_Symbol,value.var = 'value')
    rownames(x)                      <- x$tcga.id
    x$tcga.id                        <- NULL
    mutation.data                    <- as.matrix(x)
    mutation.data[mutation.data > 1] <- 1
    mutation.data
}
names(TCGA.mutation.data.list) <- names(TCGA.cancer.file.meta)

save(file = 'client-side/output/organize.TCGA.mutation.data.R.output/organize.TCGA.mutation.data.RData',list=c('TCGA.mutation.data.list'))


