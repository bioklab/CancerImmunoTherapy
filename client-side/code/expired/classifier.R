IFNg.up.gene <- res$gene.id[res$padj < 0.001 & res$log2FoldChange > 0]
IFNg.dn.gene <- res$gene.id[res$padj < 0.001 & res$log2FoldChange < 0]
require(GSVA)
require(ggplot2)
load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')



CD8A <- 'ENSG00000153563'
CD8B <- 'ENSG00000172116'
GZMA <- 'ENSG00000145649'
GZMB <- 'ENSG00000100453'
PRF1 <- 'ENSG00000180644'



load('server-side/RData/Resp.and.non.resp/SRP070710_Melanoma.RData')
l
SRP070710_Melanoma_log2.fpkm.matrix <- SRP070710_Melanoma_log2.fpkm.matrix[intersect(protein.coding.gene.ensemble.id,rownames(SRP070710_Melanoma_log2.fpkm.matrix)),]
tmp <- gsva(expr = SRP070710_Melanoma_log2.fpkm.matrix,
                 gset.idx.list = list(IFNg.up.gene=IFNg.up.gene,IFNg.dn.gene=IFNg.dn.gene),
                 method = 'ssgsea',
                 ssgsea.norm=FALSE)
score <- tmp[1,] - tmp[2,]
#score <- tmp[1,] 
#score <- -1 * tmp[2,] 

#CD8.T.cell.score <- apply(SRP070710_Melanoma_log2.fpkm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)

SRP070710_Metadata <- SRP070710_Metadata[SRP070710_Metadata$previous_mapki =='Y',]
SRP070710_Metadata$score <- score[SRP070710_Metadata$Run %>% as.character]
#SRP070710_Metadata$cd8score <- CD8.T.cell.score[SRP070710_Metadata$Run %>% as.character]

ggplot(SRP070710_Metadata,aes(x=anti_pd_1_response,y=score)) + geom_boxplot() + geom_point(aes(x=anti_pd_1_response,y=score))
wilcox.test(SRP070710_Metadata$score[SRP070710_Metadata$anti_pd_1_response == 'Complete Response'],SRP070710_Metadata$score[SRP070710_Metadata$anti_pd_1_response == 'Progressive Disease'])


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
CD8.T.cell.score <- apply(SRP094781_log2.fpkm.matrix[c(CD8A,CD8B,GZMA,GZMB,PRF1),],2,mean)

ggplot(SRP094781_Metadata,aes(x=response,y=score)) + geom_boxplot() + geom_point(aes(x=response,y=score))
wilcox.test(SRP094781_Metadata$score[SRP094781_Metadata$response == 'PRCR'],SRP094781_Metadata$score[SRP094781_Metadata$response == 'PD'])


