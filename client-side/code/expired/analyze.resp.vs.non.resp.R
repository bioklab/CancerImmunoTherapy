require(plyr)
require(dplyr)
require(ggplot2)
require(GSVA)
require(foreach)
require(quantreg)
require(DESeq2)



load('server-side/RData/Resp.and.non.resp/SRP070710_Melanoma.RData')

df                      <- SRP070710_Metadata
df                      <- df[df$previous_mapki == 'N',] # treatment naive
df$is.response          <- ifelse(grepl(x=df$anti_pd_1_response,pattern='Response'),'yes','no')
df$is.response          <- factor(x = df$is.response,levels = c('no','yes'))
rownames(df)            <- df$Run

read.count.matrix       <- round(2^SRP070710_Melanoma_log2.read.count.matrix - 1)
read.count.matrix       <- read.count.matrix[,df$Run %>% as.character]
flag                    <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix)* 0.5)
read.count.matrix       <- read.count.matrix[flag,]
dds                     <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = df[colnames(read.count.matrix),],design = ~ is.response )
dds                     <- DESeq(dds)
res                     <- results(dds,contrast = c('is.response','yes','no'))
res                     <- as.data.frame(res)
res                     <- res[complete.cases(res),]
res                     <- res[order(res$pvalue),]
res1                    <- res





#########


load('~/Project/PNMA5-PD1-response/server-side//RData//Immunotherapy_Melanoma_SRP094781.RData')

df                <- SRP094781_Metadata[SRP094781_Metadata$Treatment_Status == 'Pre',]
df                <- df[df$response %in% c('PRCR','PD'),]
df$Patient_Number <- paste('Pt',df$Patient_Number,sep='') # treatment naive
tmp               <- read.table(file = '~/Project/PNMA5-PD1-response/server-side/RData//SRP094781_Ipi_naive.txt',header=FALSE)
ipi.naive.patient <- tmp$V1
df                <- df[(df$Patient_Number  %in% ipi.naive.patient) == FALSE,]
rownames(df)      <- df$Run

read.count.matrix       <- round(2^SRP094781_log2.read.count.matrix - 1)
read.count.matrix       <- read.count.matrix[,df$Run %>% as.character]
flag                    <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix)* 0.5)
read.count.matrix       <- read.count.matrix[flag,]
dds                     <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = df[colnames(read.count.matrix),],design = ~ response )
dds                     <- DESeq(dds)
res                     <- results(dds,contrast = c('response','PRCR','PD'))
res                     <- as.data.frame(res)
res                     <- res[complete.cases(res),]
res                     <- res[order(res$pvalue),]

res2 <- res



#################

res1.sig.gene <- rownames(res1)[res1$padj < 0.01]
res2.sig.gene <- rownames(res2)[res2$padj < 0.01]


intersect(res1.sig.gene,res2.sig.gene)


################















#######################
require(org.Hs.eg.db)

tmp            <- select(org.Hs.eg.db, keys=rownames(res), columns='SYMBOL', keytype="ENSEMBL")
gene.symbol.df <- ddply(tmp,.(ENSEMBL), function(x) paste(x$SYMBOL,collapse=':'))
idx <- match(rownames(res),gene.symbol.df$ENSEMBL)
res$gene.symbol <- gene.symbol.df$V1[idx]
res <- res[order(res$pvalue),]



sig.gene <- ensemble2entrez[rownames(res)[res$padj < 0.001]]
patient.signature <- intersect(sig.gene,L1000.gene.entrez.id)
patient.signature <- setdiff(patient.signature,c('2817','4313'))

save(file='client-side/output/analyze.resp.vs.non.resp.R.output/analyze.resp.vs.non.resp.RData',list=c('patient.signature'))




# ggplot(df,aes(x=is.response,y=immune.score)) + geom_boxplot()
# 
# #TGFBR3
# plot(x=df$immune.score,y=df$expr)
# plot(x=df$immune.score,y=df$residual.expr)
