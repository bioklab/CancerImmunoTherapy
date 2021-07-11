# For each cell line, compute ssGSEA scores for geneset HALLMARK_INTERFERON_GAMMA_RESPONSE
require(GSVA)
require(dplyr)
require(data.table)
load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')
load('client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData')
load('server-side/RData/lincs/CRISPR/CRISPR.expr.RData')

A375.matrix                   <- cell.line.inst.matrix[,grepl(pattern = 'A375',x = colnames(cell.line.inst.matrix))]
cell.line.ssgsea.score.matrix <- gsva(expr = A375.matrix,
                                      gset.idx.list = list(IFNg.up.gene=IFNg.up.gene,IFNg.dn.gene=IFNg.dn.gene),
                                      method = 'ssgsea',
                                      ssgsea.norm=FALSE)


  
x <- cell.line.ssgsea.score.matrix['IFNg.up.gene',] - cell.line.ssgsea.score.matrix['IFNg.dn.gene',]
x <- sort(x,decreasing = TRUE)
inst.info <- fread(input='server-side/RData/lincs/CRISPR/instinfo.txt') %>% as.data.frame
rownames(inst.info) <- inst.info$distil_id
inst.info <- inst.info[names(x),]
View(inst.info[names(x),])

require(plyr)
inst.info$IFNg.score <- x[inst.info$distil_id %>% as.character]
batch.score <- ddply(inst.info,.(det_plate),function(x) median(x$IFNg.score[x$pert_mfc_desc == '-666'])    )
rownames(batch.score) <- batch.score$det_plate

adjusted.IFNg.score <- inst.info$IFNg.score - batch.score[inst.info$det_plate %>% as.character,'V1']
inst.info$adjusted.IFNg.score <- adjusted.IFNg.score

rs <- ddply(inst.info,.(pert_mfc_desc), function(x) data.frame(V1=median(x$adjusted.IFNg.score),mad=mad(x$adjusted.IFNg.score)))
rs <- rs[order(rs$V1,decreasing = TRUE),]
View(rs)


#######
load('server-side/RData/ENCODE/ENCODE.shRNA.RNASeq.RData')
HEPG2.log2.fpkm.matrix <- HEPG2.log2.fpkm.matrix[intersect(rownames(HEPG2.log2.fpkm.matrix),names(ensemble2entrez)),]

rownames(HEPG2.log2.fpkm.matrix) <- ensemble2entrez[rownames(HEPG2.log2.fpkm.matrix) %>% as.character]
HEPG2.shRNA.ssgsea.score.matrix <- gsva(expr = HEPG2.log2.fpkm.matrix[rownames(A375.matrix),],
                                      gset.idx.list = list(IFNg.up.gene=IFNg.up.gene,IFNg.dn.gene=IFNg.dn.gene),
                                      method = 'ssgsea',
                                      ssgsea.norm=FALSE)
x <- HEPG2.shRNA.ssgsea.score.matrix['IFNg.up.gene',] - HEPG2.shRNA.ssgsea.score.matrix['IFNg.dn.gene',]

metdata <- fread(input='server-side/RData/ENCODE/metadata.tsv') %>% as.data.frame
rownames(metdata) <- metdata[,1]
metdata <- metdata[colnames(HEPG2.log2.fpkm.matrix),]
pca.rs <- prcomp(HEPG2.log2.fpkm.matrix %>% t)
#df <- data.frame(pca.rs$x[,1:3],batch=metdata[,'Experiment date released'])

metdata$score <-x[metdata[,1] %>% as.character]
metdata$pc1 <- pca.rs$x[,1]
metdata$pc2 <- pca.rs$x[,2]

m <- metdata[,c('pc1','pc2')]

dist.matrix <- dist(m) %>% as.matrix

control.id <- metdata[grepl(x=metdata[,'Biosample genetic modifications targets'],pattern='control-human'),1] %>% as.character

dist.matrix <- dist.matrix[,control.id]

adjusted.score <- foreach(id = rownames(dist.matrix),.combine='c') %do% {
   vec <- dist.matrix[id,] %>% sort  
   c.id <- names(vec)[1]
   x[id] - x[c.id]
  
}
names(adjusted.score) <- rownames(dist.matrix)
adjusted.score <- sort(adjusted.score,decreasing = TRUE)
View(metdata[names(adjusted.score),])





require(LICORS)
cluster.rs <- kmeanspp(data=as.matrix(df[,1:3]),k = 3)
df$cluster <-cluster.rs$cluster
ggplot(df,aes(x=PC1,y=PC2,color=factor(cluster))) + geom_point()

View(metdata[names(x),])



