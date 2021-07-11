cell.line <- 'A375'
load(sprintf('client-side/output/analyze.L1000.data.R.output/positive//%s.analyze.L1000.RData',cell.line))
load(sprintf('server-side/RData/lincs/%s.lincs.RData',cell.line))
source('client-side/code/ggplot.style.R')
rownames(cell.line.inst.info.df) <- cell.line.inst.info.df$inst_id
cp.df                            <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_cp' | cell.line.inst.info.df$pert_iname == 'DMSO',]

cp <- 'GSK-1059615'  
#cp <- 'CTB'
# cp <- 'dinaciclib'
#cp <- 'tozasertib'
# cp <- 'TG-101348'
# cp <- 'pyrvinium-pamoate'
# cp <- 'pitavastatin'
cp <- 'NVP-TAE226'
cp <- 'AZD-6482'
draw.df <- cp.df[cp.df$pert_iname %in% c(cp,'DMSO'),]
draw.df[draw.df$pert_iname == 'DMSO','pert_dose'] <- 0
draw.df$score <- adjusted.IFNg.ssgsea.score[draw.df$inst_id %>% as.character]
ggplot(draw.df,aes(x=factor(pert_dose),y=score,color=factor(pert_time))) +geom_boxplot() + geom_point() +facet_wrap( ~ pert_time, nrow = 1) + stat_smooth(formula = y ~ x, method='lm') + ylab('ssGSEA.score')

ggplot(draw.df,aes(x=(pert_dose),y=score))  + geom_point()  + stat_smooth(formula = y ~ x + 0, method='lm') + ggplot.style +facet_wrap( ~ pert_time, nrow = 1)+ ylab('ssGSEA.score') + xlab('dosage (um)')

ggplot(draw.df,aes(x=(pert_dose),y=score))  + geom_point()  + stat_smooth(formula = y ~ x + 0, method='lm') + ggplot.style 



#gene <- 'SCP2'
gene <- c('lacZ','JAK2','CCND1','DUSP28','UnTrt')

oe.df <- cell.line.inst.info.df[cell.line.inst.info.df$rna_plate %in% trt.oe.plate &  cell.line.inst.info.df$pert_time == 96,]

draw.df <- oe.df[oe.df$pert_iname %in% c(gene,'UnTrt'),]
draw.df$score <- adjusted.IFNg.ssgsea.score[draw.df$inst_id %>% as.character]
ggplot(draw.df,aes(x=factor(pert_iname),y=score)) +geom_boxplot() + geom_point() + ggplot.style + xlab('') + ylab('ssGSEA.score')


# df <- read.csv(file='server-side/RData/lincs/LINCS.data.statistics.csv',header=TRUE)
# cp.cell.line <- df$cell_id[df$trt_cp > 6000] %>% as.character
# oe.cell.line <- df$cell_id[df$trt_oe > 0]    %>% as.character
# 
# cell.line.vec <- cp.cell.line
# 
# pool.rs <- foreach(cell.line = cell.line.vec,.combine='rbind') %do% {
#   
#   lincs.data.file   <- sprintf('server-side/RData/lincs/%s.lincs.RData',cell.line)
#   load(lincs.data.file)
#   cp.df <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_cp' | cell.line.inst.info.df$pert_iname == 'DMSO',]
#   rs <- ddply(cp.df,.(pert_iname,pert_time,pert_dose),nrow)
#   rs$cell.line <- cell.line
#   rs <- rs[order(rs$V1,decreasing = TRUE),]
# }
# 
# 
# 
# estimate.sigma <- function(cp.df,pert.time){
#   if(sum(cp.df$pert_time == pert.time) == 0) {return(-1)}
#   tmp <- cp.df[cp.df$pert_time == pert.time ,]     
#   rs <- ddply(tmp,.(pert_iname,pert_time,pert_dose),nrow)  
#   rs <- rs[order(rs$V1,decreasing = TRUE),]
#   pert.iname <- rs$pert_iname[1]
#   pert.dose <- rs$pert_dose[1]
#   inst.id <- cp.df$inst_id[cp.df$pert_time == pert.time & cp.df$pert_iname == pert.iname & cp.df$pert_dose == pert.dose]
#   print(sprintf('%sh %s',pert.time,pert.iname))
#   #1.4826 * mad(adjusted.IFNg.ssgsea.score[inst.id])
# }


##########
cell.line <- 'A549'

load('client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData')
load('client-side/output/analyze.PTPN2.RNASeq.data.R.output/analyze.PTPN2.RNASeq.data.RData')

input.file  <- sprintf('server-side/RData/lincs/%s.lincs.RData',cell.line)
load(input.file)

cp.name   <- 'dinaciclib'
cp.name <- 'mycophenolic-acid'
cp.name <- 'BGT-226'
ctrl.name <- 'DMSO'
pert_type <- 'trt_cp'

rna.plate.vec <- cell.line.inst.info.df$rna_plate[cell.line.inst.info.df$pert_iname == cp.name & cell.line.inst.info.df$pert_type == pert_type] %>% unique %>% as.character
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
  
  cp.inst.id       <- cell.line.inst.info.df$inst_id[cell.line.inst.info.df$rna_plate == rna.plate & cell.line.inst.info.df$pert_iname == cp.name ] %>% as.character
  cp.profile       <- cell.line.inst.matrix[,cp.inst.id]
  if(length(cp.inst.id) == 1){
      cp.profile    <- rank(cp.profile) 
      delta.profile <- cp.profile - ave.ctrl.profile 
      delta.profile <- matrix(delta.profile,ncol=1)
  }else{
      cp.profile       <- apply(cp.profile,2,rank)
      delta.profile    <- apply(cp.profile,2,function(x) x - ave.ctrl.profile)
  }

  colnames(delta.profile) <- cp.inst.id
  rownames(delta.profile) <- rownames(cell.line.inst.matrix)
  delta.profile
}
rownames(cell.line.inst.info.df) <- cell.line.inst.info.df$inst_id
pert.hour <- cell.line.inst.info.df[colnames(delta.profile),'pert_time']
pert.dose <- cell.line.inst.info.df[colnames(delta.profile),'pert_dose']




mean.profile <- apply(delta.profile[,pert.hour == 24 & pert.dose == 10.00],1,median)
mean.profile[IFNg.up.gene] %>% sort
mean.profile[IFNg.dn.gene] %>% sort



# signature.gene <- c(IFN)
# signature.gene <- as.character(signature.gene)
# 
# require(org.Hs.eg.db)
# 
# signature.gene.name <- select(org.Hs.eg.db,keys=signature.gene,keytype='ENTREZID',columns='SYMBOL')
# 
# rownames(cell.line.inst.info.df) <- cell.line.inst.info.df$inst_id
# 
# 
# 
# pdf('~/Desktop/hehe.pdf')
# for(i in 1:length(signature.gene)){
#   cp.pert.id <- cell.line.inst.info.df$inst_id[cell.line.inst.info.df$pert_iname == cp.name & cell.line.inst.info.df$pert_type == 'trt_cp' & cell.line.inst.info.df$pert_dose > 0 & cell.line.inst.info.df$pert_time >0 ]
#   df1           <- data.frame(expr=cell.line.inst.matrix[signature.gene[i],cp.pert.id],treat = cp.name,plate.id = cell.line.inst.info.df[cp.pert.id,'rna_plate'])
#   ctrl.pert.id <- cell.line.inst.info.df$inst_id[cell.line.inst.info.df$pert_iname == ctrl.name ]
#   df2           <- data.frame(expr=cell.line.inst.matrix[signature.gene[i],ctrl.pert.id],treat = ctrl.name,plate.id = cell.line.inst.info.df[ctrl.pert.id,'rna_plate'])
#   c.plate.id <- intersect(df1$plate.id,df2$plate.id) %>% unique %>% as.character
#   df2           <- df2[df2$plate.id %in% c.plate.id,]
#   df1           <- df1[df1$plate.id %in% c.plate.id,]
#   
#   p <- ggplot(rbind(df1,df2),aes(y=expr,x=treat)) + geom_boxplot() + geom_point() + facet_wrap(~plate.id) + ggtitle(signature.gene.name$SYMBOL[i])
#   print(p)
# }
# dev.off()


###########

