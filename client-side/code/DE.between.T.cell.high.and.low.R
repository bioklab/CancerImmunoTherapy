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
tumor.purity.df           <- read.csv(file = 'client-side/external.data/expired/TCGA.tumor.purity.csv',header=TRUE)
tumor.purity.df$Sample.ID <- gsub(x = tumor.purity.df$Sample.ID,pattern = '01A$',replacement = '01')
rownames(tumor.purity.df) <- tumor.purity.df$Sample.ID

CD8A <- 'ENSG00000153563'
CD8B <- 'ENSG00000172116'
GZMA <- 'ENSG00000145649'
GZMB <- 'ENSG00000100453'
PRF1 <- 'ENSG00000180644'


pan.cancer.rs <- foreach(cancer.type = names(TCGA.cancer.file.meta)) %do% {
    ################### load rna-seq data ##################################################################
    rna.seq.file <- sprintf('server-side//RData/TCGA/%s',TCGA.cancer.file.meta[[cancer.type]])  
    load(rna.seq.file)
    
    ################### compute adjusted CD8 T cell inflitration score######################################
    CD8.T.cell.score   <- apply(log2.fpkm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)
    #CD8.T.cell.score   <- 2^apply(log2.tpm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)
    
    meta.df            <- data.frame(tumor.purity=tumor.purity.df[names(CD8.T.cell.score),'CPE'],CD8.T.cell.score=CD8.T.cell.score)
    rownames(meta.df)  <- names(CD8.T.cell.score)
    meta.df            <- meta.df[complete.cases(meta.df),]
    loess.fit          <- loess(data = meta.df,formula=CD8.T.cell.score~tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
    sd.value                          <- 1.4826 * mad(loess.fit$residuals)
    md.value                          <- median(loess.fit$residuals)
    meta.df$adjusted.CD8.T.cell.score <- (loess.fit$residuals - md.value)/sd.value
    
    
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
  
    
    ################### gsea analysis, how mutation affects T cell inflitration##########
    adjusted.CD8.T.cell.score        <- meta.df$adjusted.CD8.T.cell.score
    names(adjusted.CD8.T.cell.score) <- rownames(meta.df)
    common.sample                    <- intersect(rownames(mutation.data),names(adjusted.CD8.T.cell.score))
    
    score.for.gsea                   <- adjusted.CD8.T.cell.score[common.sample]
    tmp                              <- mutation.data[common.sample,]
    s                                <- apply(tmp,2, function(x) sum(x))
    tmp                              <- tmp[,s >=2]
    signature.list <- foreach(g= colnames(tmp)) %do% {
        v <- tmp[,g]  
        rownames(tmp)[v==1] 
    }
    names(signature.list) <- colnames(tmp)
    gsea.rs               <- fgsea(stats= score.for.gsea,pathways=signature.list,minSize=2, maxSize=500, nperm=10000)
    gsea.rs               <- gsea.rs[order(gsea.rs$pval),]
    gsea.rs$cancer.type   <- cancer.type
    mutation.gsea.rs      <- gsea.rs

    ############ analyze gene expression - T cell inflitration association ##########
    median.expr    <- apply(log2.fpkm.matrix,1,median)
    expressed.gene <- names(median.expr)[median.expr > 1]
    data.for.fit   <- meta.df
    expr.matrix    <- log2.fpkm.matrix[expressed.gene,rownames(data.for.fit)]
    pval.vec <- foreach(g = rownames(expr.matrix),.combine='c') %do% {
        data.for.fit$expr <- expr.matrix[g,]  
        loess.fit         <- loess(data = data.for.fit,formula=expr~tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
        adjusted.expr     <- loess.fit$residuals
        lm.df             <- data.frame(x=adjusted.expr %>% rank,y=data.for.fit$adjusted.CD8.T.cell.score %>% rank)
        lm.fit.rs         <- lm(data = lm.df,formula =  y ~ x)
        summary(lm.fit.rs)$coefficients['x',4]
    }
    names(pval.vec) <- rownames(expr.matrix)
    rna.seq.rs <- data.frame(pval=pval.vec,cancer.type=cancer.type)

    list(rna.seq.rs=rna.seq.rs,mutation.gsea.rs=mutation.gsea.rs)
}




tmp <- foreach(item = pan.cancer.rs,.combine='rbind') %do% {
    df <- item$rna.seq.rs
    df$gene <- rownames(df)
    df
}

get.fisher.p.value <- function(df) {
  test.statistics  <- -2 * sum(log(df$pval))
  c.p.value        <- pchisq(q=test.statistics,df = 2 * nrow(df),lower.tail = FALSE) 
  c.p.value        
}

combined.p.value.df           <- ddply(tmp,.(gene),get.fisher.p.value)





mean.cor.value <- ddply(tmp,.(gene),function(x) median(x$cor.value))
mean.cor.value <- mean.cor.value[order(mean.cor.value$V1),]




get.fisher.p.value <- function(df) {
    test.statistics  <- -2 * sum(log(df$pval))
    c.p.value        <- pchisq(q=test.statistics,df = 2 * nrow(df),lower.tail = FALSE) 
    c.p.value        
}

combined.p.value.df           <- ddply(pan.cancer.mutation.gsea.rs,.(pathway),get.fisher.p.value)
colnames(combined.p.value.df) <- c('gene','p.value')
combined.p.value.df           <- combined.p.value.df[order(combined.p.value.df$p.value),]







cmd.str <- sprintf('ls %s | grep maf',folder)
file.list <- system(cmd.str,intern = TRUE)





load('server-side/RData/TCGA//Lung Squamous Cell Carcinoma.RData')







CD8A <- 'ENSG00000153563'
CD8B <- 'ENSG00000172116'
GZMA <- 'ENSG00000145649'
GZMB <- 'ENSG00000100453'
PRF1 <- 'ENSG00000180644'
CD8.T.cell.score <- apply(log2.fpkm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)

# maybe here we can also use CD4.T.cell.score - CD8.T.cell.score

#CD8.T.cell.score <- log2.fpkm.matrix['ENSG00000111537',]

tumor.purity.df           <- read.csv(file = 'client-side/external.data/expired/TCGA.tumor.purity.csv',header=TRUE)
tumor.purity.df$Sample.ID <- gsub(x = tumor.purity.df$Sample.ID,pattern = '01A$',replacement = '01')
rownames(tumor.purity.df) <- tumor.purity.df$Sample.ID

median.expr      <- apply(log2.fpkm.matrix,1,median)
expressed.gene   <- rownames(log2.fpkm.matrix)[median.expr > 1]


meta.df           <- data.frame(tumor.purity=tumor.purity.df[names(CD8.T.cell.score),'CPE'],CD8.T.cell.score=CD8.T.cell.score)
rownames(meta.df) <- names(CD8.T.cell.score)
meta.df           <- meta.df[complete.cases(meta.df),]
loess.fit         <- loess(data = meta.df,formula=CD8.T.cell.score~tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
sd.value                          <- 1.4826 * mad(loess.fit$residuals)
md.value                          <- median(loess.fit$residuals)
meta.df$adjusted.CD8.T.cell.score <- (loess.fit$residuals - md.value)/sd.value


plot(x=meta.df$tumor.purity,y=meta.df$CD8.T.cell.score,xlab='tumor.purity',ylab='CD8.T.cell.score')
lowess(x=meta.df$tumor.purity,y=meta.df$CD8.T.cell.scor) %>% lines(lwd=5) # maybe age here should also be considered

plot(x=meta.df$tumor.purity,y=meta.df$adjusted.CD8.T.cell.score,xlab='tumor.purity',ylab='adjusted.CD8.T.cell.score')
lowess(x=meta.df$tumor.purity,y=meta.df$adjusted.CD8.T.cell.score) %>% lines(lwd=5)


# sample.id  <- rownames(df)[df$x >=0.6]
# residuals  <- loess.fit$residuals[sample.id]
# residuals  <- sort(residuals)
# 
# # sample.number <- as.integer(length(residuals) * 0.25)
# # T.cell.low     <- names(residuals)[1:sample.number]
# # T.cell.high    <- names(residuals)[(length(residuals) - sample.number+1):length(residuals)]
# # meta.df        <- data.frame(sample.id= c(T.cell.low,T.cell.high))
# # meta.df$status <- 'T.cell.low'
# # meta.df[meta.df$sample.id %in% T.cell.high,'status'] <- 'T.cell.high'
# meta.df                  <- data.frame(sample.id=sample.id)
# rownames(meta.df)        <- meta.df$sample.id
# meta.df$tumor.purity     <- tumor.purity.df[rownames(meta.df) %>% as.character,'CPE']
# meta.df$CD8.T.cell.score <- residuals[rownames(meta.df) %>% as.character]
# 
# 
# 
# boxplot(cbind(tumor.purity.df[T.cell.high,'CPE'],tumor.purity.df[T.cell.low,'CPE']))
# wilcox.test(tumor.purity.df[T.cell.high,'CPE'],tumor.purity.df[T.cell.low,'CPE'])
# boxplot(cbind(CD8.T.cell.score[T.cell.high],CD8.T.cell.score[T.cell.low]))
# 
# 
# 
# meta.df           <- meta.df[complete.cases(meta.df),]


#meta.df            <- meta.df[meta.df$tumor.purity >= 0.6 ,]
meta.df$status     <- 'T.cell.low'
meta.df[meta.df$adjusted.CD8.T.cell.score >= 1.5,'status'] <- 'T.cell.high'
#meta.df[meta.df$adjusted.CD8.T.cell.score < -1.5,'status'] <- 'T.cell.low'


folder  <- 'client-side/external.data/from.gdac.firehose/gdac.broadinstitute.org_LUSC.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0'
cmd.str <- sprintf('ls %s | grep maf',folder)
file.list <- system(cmd.str,intern = TRUE)

mutation.data.df <- foreach(file = file.list,.combine='rbind') %do% {
    tmp     <- strsplit(file,split = '\\.') %>% unlist 
    tcga.id <- tmp[1]   
    #tcga.id <- gsub(x = tcga.id,pattern = '\\-',replacement = '\\.')
    file <- sprintf("%s/%s",folder,file)
    data <- fread(input = file) %>% as.data.frame
    data <- data[data$Variant_Classification != 'Silent',]
    if(nrow(data) == 0){
        return(NULL)
    }
    data$tcga.id <- tcga.id
    data[,c('tcga.id','Hugo_Symbol')]
}
mutation.data.df$value <- 1
x <- dcast(mutation.data.df,formula = tcga.id ~ Hugo_Symbol,value.var = 'value')
rownames(x) <- x$tcga.id
x$tcga.id <- NULL
mutation.data <- as.matrix(x)
mutation.data[mutation.data > 1] <- 1


mutation.data <- mutation.data[intersect(rownames(mutation.data),rownames(meta.df)),]
meta.df <- meta.df[rownames(mutation.data),]


#require(phenoTest)

mutation.freq <- apply(mutation.data,2,function(x) sum(x)/nrow(mutation.data))
mutation.data <- mutation.data[,mutation.freq>= 2/nrow(mutation.data) ]

signature.list <- foreach(g= colnames(mutation.data)) %do% {
   v <- mutation.data[,g]  
   rownames(mutation.data)[v==1] 
}

names(signature.list) <- colnames(mutation.data)
score <- meta.df$adjusted.CD8.T.cell.score
names(score) <- rownames(meta.df)
#gsea.rs <- gsea(x= score,gsets=signature.list,logScale=FALSE, absVals=FALSE,B=100000,center = TRUE)
gsea.rs <- fgsea(stats= score,pathways=signature.list,minSize=5, maxSize=500, nperm=10000)
gsea.rs <- gsea.rs[order(gsea.rs$pval),]
View(gsea.rs)



plotEnrichment(signature.list[['TP53']],stats = score.for.gsea,ticksSize = 0.5)
plotEnrichment(signature.list[['CFAP58']],stats = score) # for lung cancer, minSize=2 in fgsea function

x      <- summary(gsea.rs)
View(x[order(x[,4]),])



p.value.vec <- foreach(g=colnames(mutation.data),.combine='c') %do% {
    meta.df$mutation <- mutation.data[,g]  
    tmp              <- meta.df$adjusted.CD8.T.cell.score[meta.df$mutation == 0] # score for wild type  samples
    miu.estimate     <- median(tmp)
    sigma.estimate   <- 1.4826 * mad(tmp)
    tmp1             <- meta.df$adjusted.CD8.T.cell.score[meta.df$mutation == 1] # score for  mutation samples
    p.value          <- pnorm(q =tmp1,mean = miu.estimate,sd=sigma.estimate,lower.tail = FALSE)
    test.statistics  <- -2 * sum(log(p.value))
    c.p.value        <- pchisq(q=test.statistics,df = 2 * length(p.value),lower.tail = FALSE) 
    c.p.value  
    # well, maybe try rank test again?
}
names(p.value.vec) <- colnames(mutation.data)


p.value.vec <- sort(p.value.vec)
View(p.value.vec)
#p.value.vec <- p.value.vec[p.value.vec != -1]
g <- 'TRPM6'
meta.df$mutation <- mutation.data[rownames(meta.df),g]
ggplot(meta.df,) + geom_boxplot(aes(x=factor(mutation),y=adjusted.CD8.T.cell.score)) + geom_point(aes(x=factor(mutation),y=adjusted.CD8.T.cell.score))
table(meta.df$mutation)







cnv.file  <- 'client-side/external.data/from.gdac.firehose/gdac.broadinstitute.org_LUAD.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0/LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt'
data <- fread(input = cnv.file) %>% as.data.frame

get.TCGA.id <- function(x) {
    tmp <- strsplit(x,split = '\\-') %>% unlist
    tmp[4] <- gsub(tmp[4],pattern = 'A',replacement = '')
    paste(tmp[1:4],collapse = '-')
}
data$sample.id <- sapply(data$Sample,get.TCGA.id)
data$cohort.id <- data$sample.id
data$probe    <- paste(data$Chromosome,data$Start,data$End,sep=":")
tmp           <- data[,c('cohort.id','probe','Segment_Mean')]

cnv.data <- dcast(data = tmp,formula = cohort.id ~ probe,value.var = 'Segment_Mean',fill = 0)
rownames(cnv.data) <- cnv.data$cohort.id
cnv.data$cohort.id <- NULL
cnv.data.matrix <- as.matrix(cnv.data)
cnv.data.matrix <- cnv.data.matrix[intersect(rownames(meta.df),rownames(cnv.data.matrix)),]
meta.df         <- meta.df[intersect(rownames(meta.df),rownames(cnv.data.matrix)),]
cnv.data.matrix <- cnv.data.matrix[rownames(meta.df),]


cor.value <- foreach(g=colnames(cnv.data.matrix)[1:1000],.combine='c') %do% {
  meta.df$cnv       <- cnv.data.matrix[,g]  
  if(sum(meta.df$cnv ==0) < (0.5 * nrow(meta.df) )){
  loess.fit         <- loess(data = meta.df,formula=cnv~tumor.purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
  cor(loess.fit$residuals,meta.df$adjusted.CD8.T.cell.score,method='spearman')
  }else{
    -1
  }
  
}
names(cor.value) <- colnames(cnv.data.matrix)


h.sample <- rownames(meta.df)[meta.df$status == 'T.cell.high']
h.sample.cnv <- cnv.data.matrix[h.sample,]



high.T.cell.mutation.data <- mutation.data[intersect(rownames(mutation.data),rownames(meta.df)[meta.df$status == 'T.cell.high']),]
low.T.cell.mutation.data  <- mutation.data[intersect(rownames(mutation.data),rownames(meta.df)[meta.df$status == 'T.cell.low']),]

f.df <- foreach(g = colnames(high.T.cell.mutation.data),.combine='rbind') %do% {
    f.high <- sum(high.T.cell.mutation.data[,g]) / nrow(high.T.cell.mutation.data)  
    f.low <- sum(low.T.cell.mutation.data[,g]) / nrow(low.T.cell.mutation.data)  
    data.frame(f.high=f.high,f.low=f.low)
  
}
f.df$gene.id <- colnames(high.T.cell.mutation.data)

plot(x=f.df$f.low,y=f.df$f.high,xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1))
f.df <- f.df[order(f.df$f.high,decreasing = TRUE),]
View(f.df)
# maybe analysis rank?




read.count.matrix  <- round(2^log2.read.count.matrix - 1)
#flag              <- apply(read.count.matrix,1,function(x) sum(x>1) >= 0.5 * nrow(meta.df))
#read.count.matrix <- read.count.matrix[flag,]
read.count.matrix <- read.count.matrix[expressed.gene,rownames(meta.df) %>% as.character]

#dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df,design= ~ adjusted.CD8.T.cell.score + tumor.purity )
#dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df,design= ~ adjusted.CD8.T.cell.score)
#dds                <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df,design= ~ tumor.purity + status)
dds                <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df,design= ~ tumor.purity + status)

dds               <- DESeq(dds)
res               <- results(dds,contrast = c('status','T.cell.high','T.cell.low'))
#res               <- results(dds,name = 'adjusted.CD8.T.cell.score')
res               <- as.data.frame(res)
res               <- res[complete.cases(res),]
res$gene.id       <- rownames(res)
res               <- res[order(res$pvalue),]
res               <- res[complete.cases(res),]
sig.up.gene       <- rownames(res)[res$padj < 0.001 & res$log2FoldChange > 0]
sig.dn.gene       <- rownames(res)[res$padj < 0.001 & res$log2FoldChange < 0]
#sig.up.gene       <- rownames(res)[res$padj < 0.01 & res$log2FoldChange > 0]
#sig.dn.gene       <- rownames(res)[res$padj < 0.01 & res$log2FoldChange < 0]


########################### Trash Code ######################################
#flag                       <- rownames(log2.fpkm.matrix) %in% names(ensemble2entrez)
#log2.fpkm.matrix           <- log2.fpkm.matrix[flag,]
#rownames(log2.fpkm.matrix) <- ensemble2entrez[rownames(log2.fpkm.matrix)]
# up.res <- res[sig.up.gene,]
# up.res <- up.res[order(up.res$log2FoldChange,decreasing = TRUE),]
# 
# dn.res <- res[sig.dn.gene,]
# dn.res <- dn.res[order(dn.res$log2FoldChange,decreasing = FALSE),]
# 
# boxplot(cbind(log2.fpkm.matrix['ENSG00000169248',T.cell.high],log2.fpkm.matrix['ENSG00000169248',T.cell.low]))
# # rs.matrix <- foreach(target.gene = rownames(log2.fpkm.matrix),.combine='rbind') %do% {
# #     target.gene.expr.vec <- log2.fpkm.matrix[target.gene,]
# #     purity.vec <- tumor.purity.df[names(target.gene.expr.vec),'CPE']
# #     data.for.lm <- data.frame(score=CD8.T.cell.score[names(target.gene.expr.vec)],purity = purity.vec ,expr = target.gene.expr.vec)
# #     data.for.lm <- data.for.lm[complete.cases(data.for.lm) & data.for.lm$purity >= 0.8,]
# #     lm.rs <- lm(data.for.lm,formula = score ~ purity + expr)
# #     summary(lm.rs)$coefficients['expr',c(1,4)]
# # }
# # rownames(rs.matrix) <- rownames(log2.fpkm.matrix)
# # rs.matrix <- cbind(rs.matrix,p.adjust(rs.matrix[,2],method='fdr'))
# # colnames(rs.matrix) <- c('effect.size','p.value','fdr')
# # flag <- rs.matrix[,3] < 0.01
# # plot(x=rs.matrix[,'effect.size'],y=-1 * rs.matrix[,'fdr'] %>% log10)
# # 
# # sig.matrix <- rs.matrix[flag,]
# 
# # Below is specific for prostate cancer
# #spe.sample: look at CD.8.T.cell.score vs purity plot, the samples on the left are spe.sample
# flag <- ifelse(colnames(log2.fpkm.matrix) %in% spe.sample,'1','0')
# df   <- data.frame(f=flag,purity=tumor.purity.df[colnames(log2.fpkm.matrix),'CPE']) 
# rownames(df) <- colnames(log2.fpkm.matrix)
# df   <- df[complete.cases(df),]
# df <- df[df$purity <= 0.7,]
# 
# read.count.matrix <- round(2^log2.read.count.matrix[,rownames(df) %>% as.character] - 1)
# #flag              <- apply(read.count.matrix,1,function(x) sum(x>1) >= 0.5 * nrow(meta.df))
# #read.count.matrix <- read.count.matrix[flag,]
# read.count.matrix <- read.count.matrix[expressed.gene,]
# read.count.matrix <- read.count.matrix[expressed.gene,rownames(df)]
# 
# dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = df,design= ~ f + purity )
# dds               <- DESeq(dds)
# res               <- results(dds,contrast = c('f','1','0'))
# res               <- as.data.frame(res)
# res               <- res[complete.cases(res),]
# res$gene.id       <- rownames(res)
# res               <- res[order(res$pvalue),]
# res               <- res[complete.cases(res),]
# up.res <- res[sig.up.gene,]
# up.res <- up.res[order(up.res$log2FoldChange,decreasing = TRUE),]
# sig.up.gene       <- rownames(res)[res$padj < 0.001 & res$log2FoldChange > 0]
# sig.dn.gene       <- rownames(res)[res$padj < 0.001 & res$log2FoldChange < 0]
# 
# dn.res <- res[sig.dn.gene,]
# dn.res <- dn.res[order(dn.res$log2FoldChange,decreasing = FALSE),]

