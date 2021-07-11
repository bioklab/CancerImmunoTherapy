require(data.table)
read.count <- fread(input='client-side/external.data/GEO.dataset.meta/GSE109933_yfp_Combined.Allcounts.txt') %>% as.data.frame
#read.count <- fread(input='client-side/external.data/GEO.dataset.meta/GSE109971_bulk_Combined.Allcounts.txt') %>% as.data.frame

rownames(read.count) <- read.count$id
read.count$symbol <- NULL
read.count$id <- NULL
read.count.matrix <- as.matrix(read.count)
sample.source     <- sapply(colnames(read.count.matrix),function(x) strsplit(x=x,split = '_') %>% unlist)[2,]
meta.df           <- data.frame(sample.id = colnames(read.count.matrix),sample.source=sample.source)
rownames(meta.df) <- meta.df$sample.id
T.cell.high       <- c('7160c2','6620c1','6499c4','6555c3','6556c4','6421c2','2838c3')
T.cell.low        <- c('6556c3','6499c3','6694c2','6422c1','2699c4','6419c1','6419c5','5821c3')
read.count.matrix <- read.count.matrix[,sample.source %in% c(T.cell.high,T.cell.low)]

meta.df$status   <- 'haha'
meta.df[meta.df$sample.source %in% T.cell.high,'status'] <- 'T.cell.high'
meta.df[meta.df$sample.source %in% T.cell.low,'status'] <- 'T.cell.low'
meta.df <- meta.df[meta.df$status != 'haha',]

flag <- apply(read.count.matrix,1,function(x) sum(x > 1) >= 5)
read.count.matrix <- read.count.matrix[flag,]

require(DESeq2)
dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df[colnames(read.count.matrix),],design= ~ status)
dds               <- DESeq(dds)
res               <- results(dds,contrast = c('status','T.cell.high','T.cell.low'))
res               <- as.data.frame(res)
res               <- res[complete.cases(res),]
res$gene.id       <- rownames(res)
res               <- res[order(res$padj),]
View(res[res$padj < 0.05 & res$log2FoldChange < 0,])
View(res[res$padj < 0.05 & res$log2FoldChange > 0,])
