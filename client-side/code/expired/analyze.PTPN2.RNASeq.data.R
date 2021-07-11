# Analyze PTPN2 CRISPR-RNASeq data. DE between IFNg vs NS, PTPN2 sgRNA
require(DESeq2)
require(data.table)
require(plyr)
require(dplyr)
require(foreach)
load('server-side/RData/PTPN2.sgRNA/PTPN2.sgRNA.expression.RData')
load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')

IFN.gamma.response.gene             <- intersect(L1000.gene.entrez.id,genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']])
idx                                 <- match(IFN.gamma.response.gene,ensemble2entrez)
IFN.gamma.response.gene.ensemble.id <- names(ensemble2entrez)[idx]

idx                                 <- match(L1000.gene.entrez.id,ensemble2entrez)
L1000.gene.ensemble.id              <- names(ensemble2entrez)[idx]



meta.df                    <- fread(input='client-side/external.data/GEO.dataset.meta/GSE99299.exp.meta.txt',header=TRUE) %>% as.data.frame
meta.df[meta.df$cell_line == 'HT29','cell_type'] <- 'colon cancer'
rownames(meta.df)          <- meta.df$Run
meta.df$genotype_variation <- relevel(x=meta.df$genotype_variation %>% factor,ref = 'Control')
meta.df                    <- meta.df[meta.df$Organism != 'Mus musculus' &  meta.df$genotype_variation == 'Ptpn2ko'   ,]

#meta.df                    <- meta.df[meta.df$genotype_variation != 'NS',]
#meta.df                    <- meta.df[meta.df$treated_with == 'IFNg' & meta.df$cell_line == 'MelJuSo',]

# pdf('client-side/output/analyze.PTPN2.RNASeq.data.R.output/IFNg.expr.de.pdf')
# for(gene in IFN.gamma.response.gene.ensemble.id){
#     meta.df$expr <- log2.fpkm.matrix[gene,rownames(meta.df)]  
#     p <- ggplot(meta.df,aes(x=treated_with,y=expr,col=treated_with)) + geom_point(size=5) + facet_wrap(~cell_line) + ggtitle(gene)
#     print(p)
# }
# dev.off()


srr.id            <- meta.df$Run %>% as.character
read.count.matrix <- round(2^log2.read.count.matrix[,srr.id] -1) 
flag              <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix) * 0.5)
read.count.matrix <- read.count.matrix[flag,]             
dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df[colnames(read.count.matrix),],design= ~ treated_with)
dds               <- DESeq(dds)
res               <- results(dds,contrast = c('treated_with','IFNg','NS'))
res               <- as.data.frame(res)
res               <- res[complete.cases(res),]
res$gene.id       <- rownames(res)

rs <- res[L1000.gene.ensemble.id,]
rs <- rs[order(rs$pvalue),]
rs <- rs[complete.cases(rs),]

IFNg.up.gene <- rs$gene.id[rs$padj < 0.001 & rs$log2FoldChange > 0]
IFNg.dn.gene <- rs$gene.id[rs$padj < 0.001 & rs$log2FoldChange < 0]
IFNg.up.gene <- ensemble2entrez[IFNg.up.gene] %>% as.character
IFNg.dn.gene <- ensemble2entrez[IFNg.dn.gene] %>% as.character
save(file='client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData',list=c('IFNg.up.gene','IFNg.dn.gene'))

ptpn2.res          <- res
IFNg.up.gene.ptpn2 <- res$gene.id[res$padj < 0.001 & res$log2FoldChange > 0]
IFNg.dn.gene.ptpn2 <- res$gene.id[res$padj < 0.001 & res$log2FoldChange < 0]

#IFNg.up.gene <- rs$gene.id[rs$padj < 0.001 & rs$log2FoldChange > 0]
#IFNg.dn.gene <- rs$gene.id[rs$padj < 0.001 & rs$log2FoldChange < 0]
pca.rs <- prcomp(log2.fpkm.matrix[c(IFNg.up.gene.ptpn2,IFNg.dn.gene.ptpn2),rownames(meta.df)] %>% t)
plot(pca.rs$x[,1:2])
meta.df$pc1 <- pca.rs$x[,1]
meta.df$pc2 <- pca.rs$x[,2]
ggplot(meta.df,aes(x=pc1,y=pc2,col=treated_with,shape=cell_line)) + geom_point(size=5)



############ IFNg vs non-treatment, wildtype
meta.df                    <- fread(input='client-side/external.data/GEO.dataset.meta/GSE99299.exp.meta.txt',header=TRUE) %>% as.data.frame
meta.df[meta.df$cell_line == 'HT29','cell_type'] <- 'colon cancer'
rownames(meta.df)          <- meta.df$Run
meta.df$genotype_variation <- relevel(x=meta.df$genotype_variation %>% factor,ref = 'Control')
meta.df                    <- meta.df[meta.df$Organism != 'Mus musculus' &  meta.df$genotype_variation != 'Ptpn2ko'   ,]



srr.id            <- meta.df$Run %>% as.character
read.count.matrix <- round(2^log2.read.count.matrix[,srr.id] -1) 
flag              <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix) * 0.5)
read.count.matrix <- read.count.matrix[flag,]             
dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df[colnames(read.count.matrix),],design= ~ treated_with)
dds               <- DESeq(dds)
res               <- results(dds,contrast = c('treated_with','IFNg','NS'))
res               <- as.data.frame(res)
res               <- res[complete.cases(res),]
res$gene.id       <- rownames(res)
wildtype.res      <- res


IFNg.up.gene.wildtype <- res$gene.id[res$padj < 0.001 & res$log2FoldChange > 0]
IFNg.dn.gene.wildtype <- res$gene.id[res$padj < 0.001 & res$log2FoldChange < 0]


###########################
plot(x= -1 * log10(wildtype.res[IFNg.dn.gene.ptpn2,'padj']),-1 * log10(ptpn2.res[IFNg.dn.gene.ptpn2,'padj']),xlim=c(0,20),ylim=c(0,20))
plot(x= -1 * log10(wildtype.res[IFNg.up.gene.ptpn2,'padj']),-1 * log10(ptpn2.res[IFNg.up.gene.ptpn2,'padj']),xlim=c(0,20),ylim=c(0,20))

#################
intersect(IFNg.dn.gene.ptpn2,IFNg.dn.gene.wildtype)
intersect(IFNg.up.gene.ptpn2,IFNg.up.gene.wildtype)



##################################################################################################################################################################
require(GSVA)
load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')
load('server-side/RData/Resp.and.non.resp/SRP070710_Melanoma.RData')

SRP070710_Melanoma_log2.fpkm.matrix <- SRP070710_Melanoma_log2.fpkm.matrix[intersect(protein.coding.gene.ensemble.id,rownames(SRP070710_Melanoma_log2.fpkm.matrix)),]
tmp <- gsva(expr = SRP070710_Melanoma_log2.fpkm.matrix,
            gset.idx.list = list(IFNg.up.gene=IFNg.up.gene.ptpn2,IFNg.dn.gene=IFNg.dn.gene.ptpn2),
            method = 'ssgsea',
            ssgsea.norm=FALSE)
score <- tmp[1,] - tmp[2,]
#score <- tmp[1,] 
#score <- -1 * tmp[2,] 

#CD8.T.cell.score <- apply(SRP070710_Melanoma_log2.fpkm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)

#SRP070710_Metadata <- SRP070710_Metadata[SRP070710_Metadata$previous_mapki =='Y',]
SRP070710_Metadata$score <- score[SRP070710_Metadata$Run %>% as.character]
#SRP070710_Metadata$cd8score <- CD8.T.cell.score[SRP070710_Metadata$Run %>% as.character]
colnames(SRP070710_Metadata)[2] <- 'response'
ggplot(SRP070710_Metadata,aes(x=response,y=score)) + geom_boxplot() + geom_point(aes(x=response,y=score))
wilcox.test(SRP070710_Metadata$score[SRP070710_Metadata$response == 'Complete Response'],SRP070710_Metadata$score[SRP070710_Metadata$response == 'Progressive Disease'])


load('~/Project/PNMA5-PD1-response/server-side//RData//Immunotherapy_Melanoma_SRP094781.RData')
SRP094781_log2.fpkm.matrix <- SRP094781_log2.fpkm.matrix[intersect(protein.coding.gene.ensemble.id,rownames(SRP094781_log2.fpkm.matrix)),]
SRP094781_Metadata$Patient_Number <- paste('Pt',SRP094781_Metadata$Patient_Number,sep='') # treatment naive

SRP094781_Metadata                <- SRP094781_Metadata[SRP094781_Metadata$Treatment_Status == 'Pre',]
tmp               <- read.table(file = '~/Project/PNMA5-PD1-response/server-side/RData//SRP094781_Ipi_naive.txt',header=FALSE)
ipi.naive.patient <- tmp$V1 %>% as.character
SRP094781_Metadata                <- SRP094781_Metadata[(SRP094781_Metadata$Patient_Number  %in% ipi.naive.patient) == FALSE,]




tmp <- gsva(expr = SRP094781_log2.fpkm.matrix,
            gset.idx.list = list(IFNg.up.gene=IFNg.up.gene,IFNg.dn.gene=IFNg.dn.gene),
            method = 'ssgsea',
            ssgsea.norm=FALSE)
score <- tmp[1,] - tmp[2,]


SRP094781_Metadata$score <- score[SRP094781_Metadata$Run %>% as.character]

ggplot(SRP094781_Metadata,aes(x=response,y=score)) + geom_boxplot() + geom_point(aes(x=response,y=score))
wilcox.test(SRP094781_Metadata$score[SRP094781_Metadata$response == 'PRCR'],SRP094781_Metadata$score[SRP094781_Metadata$response == 'PD'])














# IFNg.up.gene <- res$gene.id[res$padj < 0.001 & res$log2FoldChange > 0]
# IFNg.dn.gene <- res$gene.id[res$padj < 0.01 & res$log2FoldChange < 0]

# #pca.rs <- prcomp(log2.fpkm.matrix[IFN.gamma.response.gene.ensemble.id,rownames(meta.df)[meta.df$cell_line == 'MelJuSo']] %>% t)
# pca.rs <- prcomp(log2.fpkm.matrix[IFN.gamma.response.gene.ensemble.id,rownames(meta.df)] %>% t)
# 
# plot(pca.rs$x[,1:2])
# meta.df$pc1 <- pca.rs$x[,1]
# meta.df$pc2 <- pca.rs$x[,2]
# require(ggplot2)
# ggplot(meta.df,aes(x=pc1,y=pc2,col=treated_with,shape=cell_line)) + geom_point(size=5)
# 
# 
# 
# cell.line.vec <- meta.df$cell_line %>% as.character %>% unique
# rs <- foreach(cell.line = cell.line.vec,.combine='rbind') %do% {
#     flag              <- meta.df$cell_line == cell.line
#     srr.id            <- meta.df$Run[flag] %>% as.character
#     read.count.matrix <- round(2^log2.read.count.matrix[,srr.id] -1) 
#     flag              <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix) * 0.5)
#     read.count.matrix <- read.count.matrix[flag,]
#     dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df[colnames(read.count.matrix),],design= ~ genotype_variation)
#     dds               <- DESeq(dds)
#     res               <- results(dds)
#     res               <- as.data.frame(res)
#     res               <- res[complete.cases(res),]
#     res$gene.id       <- rownames(res)
#     res$cell.line     <- cell.line
#     res[IFN.gamma.response.gene.ensemble.id,]
#     res[L1000.gene.ensemble.id,]
# }
# rs <- rs[complete.cases(rs),]
# sig.dn.df <- rs[rs$padj < 0.01 & rs$log2FoldChange < 0,]
# sig.up.df <- rs[rs$padj < 0.01 & rs$log2FoldChange > 0,]
# 
# tmp.dn <- dcast(sig.dn.df,formula = gene.id ~ cell.line,value.var = 'log2FoldChange')
# 
# tmp.up <- dcast(sig.up.df,formula = gene.id ~ cell.line,value.var = 'log2FoldChange')
# rownames(tmp.up) <- tmp.up$gene.id
# tmp.up$gene.id <- NULL
# tmp.up          <- as.matrix(tmp.up)
# v <- apply(tmp.up,1,function(x) sum(is.na(x) == FALSE))
# up.gene <- names(v)[v >=2]
# PTPN2.del.up.gene <- ensemble2entrez[up.gene] %>% c %>% as.character
# 
# PTPN2.del.up.gene <- setdiff(PTPN2.del.up.gene,IFN.gamma.response.gene) %>% unique %>% c
# 
# tmp1 <- setdiff(L1000.gene.entrez.id,PTPN2.del.up.gene)
# tmp1 <- setdiff(tmp1,genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']])
# 
# random.gene <- sample(x=tmp1,size=length(PTPN2.del.up.gene))
# save(file='client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData',list=c('PTPN2.del.up.gene','random.gene'))






#dn.df <- rs[rs$log2FoldChange<0 & rs$padj < 0.01,] 










# # Analyze PTPN2 CRISPR-RNASeq data. DE between IFNg vs NS, both PTPN2 and control sgRNA
# require(DESeq2)
# require(data.table)
# require(plyr)
# require(dplyr)
# require(foreach)
# load('server-side/RData/PTPN2.sgRNA/PTPN2.sgRNA.expression.RData')
# load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')
# 
# IFN.gamma.response.gene             <- intersect(L1000.gene.entrez.id,genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']])
# idx                                 <- match(IFN.gamma.response.gene,ensemble2entrez)
# IFN.gamma.response.gene.ensemble.id <- names(ensemble2entrez)[idx]
# 
# idx                                 <- match(L1000.gene.entrez.id,ensemble2entrez)
# L1000.gene.ensemble.id              <- names(ensemble2entrez)[idx]
# 
# 
# 
# meta.df                    <- fread(input='client-side/external.data/GEO.dataset.meta/GSE99299.exp.meta.txt',header=TRUE) %>% as.data.frame
# meta.df[meta.df$cell_line == 'HT29','cell_type'] <- 'colon cancer'
# rownames(meta.df)          <- meta.df$Run
# meta.df$genotype_variation <- relevel(x=meta.df$genotype_variation %>% factor,ref = 'Control')
# meta.df                    <- meta.df[meta.df$Organism != 'Mus musculus',]
# 
# #meta.df                    <- meta.df[meta.df$genotype_variation != 'NS',]
# #meta.df                    <- meta.df[meta.df$treated_with == 'IFNg' & meta.df$cell_line == 'MelJuSo',]
# 
# # pdf('client-side/output/analyze.PTPN2.RNASeq.data.R.output/IFNg.expr.de.pdf')
# # for(gene in IFN.gamma.response.gene.ensemble.id){
# #     meta.df$expr <- log2.fpkm.matrix[gene,rownames(meta.df)]  
# #     p <- ggplot(meta.df,aes(x=treated_with,y=expr,col=treated_with)) + geom_point(size=5) + facet_wrap(~cell_line) + ggtitle(gene)
# #     print(p)
# # }
# # dev.off()
# 
# 
# srr.id            <- meta.df$Run %>% as.character
# read.count.matrix <- round(2^log2.read.count.matrix[,srr.id] -1) 
# flag              <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix) * 0.5)
# read.count.matrix <- read.count.matrix[flag,]             
# dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df[colnames(read.count.matrix),],design= ~ treated_with)
# dds               <- DESeq(dds)
# res               <- results(dds,contrast = c('treated_with','IFNg','NS'))
# res               <- as.data.frame(res)
# res               <- res[complete.cases(res),]
# res$gene.id       <- rownames(res)
# 
# rs <- res[L1000.gene.ensemble.id,]
# rs <- rs[order(rs$pvalue),]
# rs <- rs[complete.cases(rs),]
# 
# IFNg.up.gene <- rs$gene.id[rs$padj < 0.001 & rs$log2FoldChange > 0]
# IFNg.dn.gene <- rs$gene.id[rs$padj < 0.001 & rs$log2FoldChange < 0]
# IFNg.up.gene <- ensemble2entrez[IFNg.up.gene]
# IFNg.dn.gene <- ensemble2entrez[IFNg.dn.gene]
# save(file='client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData',list=c('IFNg.up.gene','IFNg.dn.gene'))
# 







# #pca.rs <- prcomp(log2.fpkm.matrix[IFN.gamma.response.gene.ensemble.id,rownames(meta.df)[meta.df$cell_line == 'MelJuSo']] %>% t)
# pca.rs <- prcomp(log2.fpkm.matrix[IFN.gamma.response.gene.ensemble.id,rownames(meta.df)] %>% t)
# 
# plot(pca.rs$x[,1:2])
# meta.df$pc1 <- pca.rs$x[,1]
# meta.df$pc2 <- pca.rs$x[,2]
# require(ggplot2)
# ggplot(meta.df,aes(x=pc1,y=pc2,col=treated_with,shape=cell_line)) + geom_point(size=5)
# 
# 
# 
# cell.line.vec <- meta.df$cell_line %>% as.character %>% unique
# rs <- foreach(cell.line = cell.line.vec,.combine='rbind') %do% {
#     flag              <- meta.df$cell_line == cell.line
#     srr.id            <- meta.df$Run[flag] %>% as.character
#     read.count.matrix <- round(2^log2.read.count.matrix[,srr.id] -1) 
#     flag              <- apply(read.count.matrix,1,function(x) sum(x>=1) > ncol(read.count.matrix) * 0.5)
#     read.count.matrix <- read.count.matrix[flag,]
#     dds               <- DESeqDataSetFromMatrix(countData = read.count.matrix,colData = meta.df[colnames(read.count.matrix),],design= ~ genotype_variation)
#     dds               <- DESeq(dds)
#     res               <- results(dds)
#     res               <- as.data.frame(res)
#     res               <- res[complete.cases(res),]
#     res$gene.id       <- rownames(res)
#     res$cell.line     <- cell.line
#     res[IFN.gamma.response.gene.ensemble.id,]
#     res[L1000.gene.ensemble.id,]
# }
# rs <- rs[complete.cases(rs),]
# sig.dn.df <- rs[rs$padj < 0.01 & rs$log2FoldChange < 0,]
# sig.up.df <- rs[rs$padj < 0.01 & rs$log2FoldChange > 0,]
# 
# tmp.dn <- dcast(sig.dn.df,formula = gene.id ~ cell.line,value.var = 'log2FoldChange')
# 
# tmp.up <- dcast(sig.up.df,formula = gene.id ~ cell.line,value.var = 'log2FoldChange')
# rownames(tmp.up) <- tmp.up$gene.id
# tmp.up$gene.id <- NULL
# tmp.up          <- as.matrix(tmp.up)
# v <- apply(tmp.up,1,function(x) sum(is.na(x) == FALSE))
# up.gene <- names(v)[v >=2]
# PTPN2.del.up.gene <- ensemble2entrez[up.gene] %>% c %>% as.character
# 
# PTPN2.del.up.gene <- setdiff(PTPN2.del.up.gene,IFN.gamma.response.gene) %>% unique %>% c
# 
# tmp1 <- setdiff(L1000.gene.entrez.id,PTPN2.del.up.gene)
# tmp1 <- setdiff(tmp1,genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']])
# 
# random.gene <- sample(x=tmp1,size=length(PTPN2.del.up.gene))
# save(file='client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData',list=c('PTPN2.del.up.gene','random.gene'))






#dn.df <- rs[rs$log2FoldChange<0 & rs$padj < 0.01,] 


