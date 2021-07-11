require(plyr)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(foreach)
require(reshape2)
require(data.table)
require(stringr)
load('client-side/output/analyze.L1000.data.R.output/cell.line.analyzed.RData')


cut.off <- 0.1
positive.oe.exp.list <- foreach(cell.line = oe.cell.line) %do% {
  file         <- sprintf('client-side/output/analyze.L1000.data.R.output/positive/%s.analyze.L1000.RData',cell.line)
  load(file)
  oe.fit.df$fdr <- p.adjust(oe.fit.df$p.value,method='fdr')
  sig.df <- oe.fit.df[oe.fit.df$fdr < cut.off & oe.fit.df$effect.size > 0,]
  sig.df <- sig.df[order(sig.df$effect.size,decreasing = TRUE),]
  flag <- grepl(x=sig.df$pert_iname,pattern = 'BRD')
  sig.df <- sig.df[!flag,]
  #head(sig.df,10)
  sig.df
}
names(positive.oe.exp.list) <- oe.cell.line

df <- sapply(positive.oe.exp.list,function(x) x$pert_iname) %>% unlist %>% table %>% as.data.frame
df <- df[order(df$Freq,decreasing = TRUE),]
colnames(df) <- c('pert_iname','Freq')

tmp <- foreach(i = names(positive.oe.exp.list),.combine='rbind') %do% {
    x  <- positive.oe.exp.list[[i]]
    x$cell.line <- i
    x
}

ave.effect.size <- ddply(tmp,.(pert_iname),function(x) mean(x$effect.size))
ave.effect.size <- ave.effect.size[order(ave.effect.size$V1,decreasing = TRUE),]

merge.df <- merge(ave.effect.size,df)
merge.df <- merge.df[order(merge.df$V1,decreasing = TRUE),]





#### Pick out genes to run qPCR, OE perturbagen######
####################################################################################################################################################
library(org.Hs.eg.db)
library(annotate)
cell.line <- 'HEPG2'

load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')
load('client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData')

input.file  <- sprintf('server-side/RData/lincs/%s.lincs.RData',cell.line)
load(input.file)

#cp.name   <- 'dinaciclib'
# cp.name <- 'mycophenolic-acid'
# cp.name <- 'entinostat'
# cp.name <- 'nilotinib'
# oe.name <- 'DUSP28'
# oe.name <- 'UFM1'
# oe.name <- 'MAGEB6'
oe.name <- 'PSMD5'

ctrl.name <- 'UnTrt'
pert_type <- 'trt_oe'

rna.plate.vec <- cell.line.inst.info.df$rna_plate[cell.line.inst.info.df$pert_iname == oe.name & cell.line.inst.info.df$pert_type == pert_type] %>% unique %>% as.character
delta.profile <- foreach(rna.plate = rna.plate.vec,.combine='cbind') %do% {
  ctrl.inst.id     <- cell.line.inst.info.df$inst_id[cell.line.inst.info.df$rna_plate == rna.plate & cell.line.inst.info.df$pert_iname == ctrl.name ] %>% as.character
  ctrl.profile     <- cell.line.inst.matrix[,ctrl.inst.id]
  if(length(ctrl.inst.id) > 1){
    ctrl.profile     <- apply(ctrl.profile,2,rank)
    ave.ctrl.profile <- apply(ctrl.profile,1,mean) %>% c
  }else{
    ave.ctrl.profile           <- rank(ctrl.profile) %>% c
    ave.ctrl.profile           <- ctrl.profile %>% c
  }
  
  oe.inst.id       <- cell.line.inst.info.df$inst_id[cell.line.inst.info.df$rna_plate == rna.plate & cell.line.inst.info.df$pert_iname == oe.name ] %>% as.character
  oe.profile       <- cell.line.inst.matrix[,oe.inst.id]
  if(length(oe.inst.id) == 1){
    oe.profile    <- rank(oe.profile) 
    delta.profile <- oe.profile - ave.ctrl.profile 
    delta.profile <- matrix(delta.profile,ncol=1)
  }else{
    oe.profile       <- apply(oe.profile,2,rank)
    delta.profile    <- apply(oe.profile,2,function(x) x - ave.ctrl.profile)
  }
  
  colnames(delta.profile) <- oe.inst.id
  rownames(delta.profile) <- rownames(cell.line.inst.matrix)
  delta.profile
}
rownames(cell.line.inst.info.df) <- cell.line.inst.info.df$inst_id
pert.hour <- cell.line.inst.info.df[colnames(delta.profile),'pert_time']
pert.dose <- cell.line.inst.info.df[colnames(delta.profile),'pert_dose']



mean.profile <- apply(delta.profile,1,median)
IFNg.up.gene.profile <- mean.profile[IFNg.up.gene] %>% sort(decreasing = TRUE)
IFNg.dn.gene.profile <- mean.profile[IFNg.dn.gene] %>% sort



names(IFNg.up.gene.profile) <- lookUp(names(IFNg.up.gene.profile), 'org.Hs.eg', 'SYMBOL') %>% unlist
names(IFNg.dn.gene.profile) <- lookUp(names(IFNg.dn.gene.profile), 'org.Hs.eg', 'SYMBOL') %>% unlist
View(IFNg.up.gene.profile)
View(IFNg.dn.gene.profile)



