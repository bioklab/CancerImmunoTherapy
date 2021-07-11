require(plyr)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(foreach)
require(reshape2)

load('client-side/output/analyze.L1000.data.R.output/cell.line.analyzed.RData')
#positive oe
pan.can.oe.fit.df <- foreach(cell.line = oe.cell.line,.combine='rbind') %do% {
  file         <- sprintf('client-side/output/analyze.L1000.data.R.output/positive/%s.analyze.L1000.RData',cell.line)
  load(file)
  oe.fit.df$fdr <- p.adjust(oe.fit.df$p.value,method='fdr')
  oe.fit.df$cell.line <- cell.line
  oe.fit.df    
  
}
pan.can.oe.fit.df$oe.name <- pan.can.oe.fit.df$pert_iname
cut.off                   <- 0.05

pan.can.oe.fit.wide.df <- dcast(pan.can.oe.fit.df,oe.name ~ cell.line,value.var='fdr')
sig.oe                 <- pan.can.oe.fit.df$oe.name[pan.can.oe.fit.df$fdr < cut.off] %>% unique %>% as.character # significant in at least one cell line
sig.oe.across.cell.line           <- pan.can.oe.fit.wide.df[pan.can.oe.fit.wide.df$oe.name %in% sig.oe,]
rownames(sig.oe.across.cell.line) <- sig.oe.across.cell.line$oe.name
sig.oe.across.cell.line$oe.name   <- NULL
sig.oe.across.cell.line           <- as.matrix(sig.oe.across.cell.line)
sig.cell.line.cnt.per.oe          <- apply(sig.oe.across.cell.line,1,function(x){
  y <- x[is.na(x) == FALSE]
  sum(y < cut.off)
})
sig.cell.line.cnt.per.oe          <- sort(sig.cell.line.cnt.per.oe,decreasing = TRUE)
positive.sig.cell.line.cnt.per.oe <- sig.cell.line.cnt.per.oe
positive.sig.oe.acorss.cell.line  <- sig.oe.across.cell.line


# ### negative oe
# pan.can.oe.fit.df <- foreach(cell.line = oe.cell.line,.combine='rbind') %do% {
#   file         <- sprintf('client-side/output/analyze.L1000.data.R.output/negative/%s.analyze.L1000.RData',cell.line)
#   load(file)
#   oe.fit.df$fdr <- p.adjust(oe.fit.df$p.value,method='fdr')
#   oe.fit.df$cell.line <- cell.line
#   oe.fit.df    
#   
# }
# pan.can.oe.fit.df$oe.name <- pan.can.oe.fit.df$pert_iname
# cut.off                   <- 0.05
# 
# pan.can.oe.fit.wide.df <- dcast(pan.can.oe.fit.df,oe.name ~ cell.line,value.var='fdr')
# sig.oe                 <- pan.can.oe.fit.df$oe.name[pan.can.oe.fit.df$fdr < cut.off] %>% unique %>% as.character # significant in at least one cell line
# sig.oe.across.cell.line           <- pan.can.oe.fit.wide.df[pan.can.oe.fit.wide.df$oe.name %in% sig.oe,]
# rownames(sig.oe.across.cell.line) <- sig.oe.across.cell.line$oe.name
# sig.oe.across.cell.line$oe.name   <- NULL
# sig.oe.across.cell.line           <- as.matrix(sig.oe.across.cell.line)
# sig.cell.line.cnt.per.oe          <- apply(sig.oe.across.cell.line,1,function(x){
#   y <- x[is.na(x) == FALSE]
#   sum(y < cut.off)
# })
# sig.cell.line.cnt.per.oe          <- sort(sig.cell.line.cnt.per.oe,decreasing = TRUE)
# negative.sig.cell.line.cnt.per.oe <- sig.cell.line.cnt.per.oe
# negative.sig.oe.acorss.cell.line  <- sig.oe.across.cell.line
# 



save(file='client-side/output/pan.cancer.oe.R.output/pan.cancer.oe.RData',list = c('positive.sig.cell.line.cnt.per.oe','positive.sig.oe.acorss.cell.line'))








################## Trash ##########################
r.idx <- apply(positive.effect.cp.across.cell.line,1,function(x) sum(x<cut.off & is.na(x) == FALSE)) %>% sort(decreasing = TRUE)
positive.effect.cp.across.cell.line <- positive.effect.cp.across.cell.line[names(r.idx),]





pan.can.cp.fit.df <- pan.can.cp.fit.df[complete.cases(pan.can.cp.fit.df),]
sig.cp.name       <- pan.can.cp.fit.df$cp.name %>% as.character %>% unique

pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
  file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
  load(file)
  sig.df <- cp.fit.df[cp.fit.df$cp.name %in% sig.cp.name,]
  if(sum(cp.fit.df$cp.name %in% sig.cp.name) == 0){
    data.frame(cp.name=NA,slope=NA,p.value=NA,p.value.negative = NA,fdr=NA,fdr.negative=NA,cell.line=cell.line)  
  }else{
    sig.df$cell.line <- cell.line
    sig.df[sig.df$fdr > cut.off,'slope'] <- 0
    sig.df
  }
}



cp.across.cell.line           <- dcast(pan.can.cp.fit.df,cell.line ~ cp.name,value.var='slope') 
rownames(cp.across.cell.line) <- cp.across.cell.line$cell.line
cp.across.cell.line$cell.line <- NULL
cp.across.cell.line           <- as.matrix(cp.across.cell.line)
tmp                           <- apply(cp.across.cell.line,2, function(x) sum( x != 0 & is.na(x) == FALSE )  )
cp.across.cell.line           <- cp.across.cell.line[,order(tmp,decreasing = TRUE)]
positive.effect.cp.across.cell.line <- cp.across.cell.line






pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
  file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
  load(file)
  sig.df <- cp.fit.df[cp.fit.df$fdr.negative < cut.off,]
  if(sum(cp.fit.df$fdr.negative < cut.off) == 0){
    data.frame(cp.name=NA,slope=NA,p.value=NA,p.value.negative = NA,fdr=NA,fdr.negative=NA,cell.line=cell.line)  
  }else{
    sig.df$cell.line <- cell.line
    sig.df
  }
}
pan.can.cp.fit.df <- pan.can.cp.fit.df[complete.cases(pan.can.cp.fit.df),]
sig.cp.name <- pan.can.cp.fit.df$cp.name %>% as.character %>% unique

pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
  file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
  load(file)
  sig.df <- cp.fit.df[cp.fit.df$cp.name %in% sig.cp.name,]
  if(sum(cp.fit.df$cp.name %in% sig.cp.name) == 0){
    data.frame(cp.name=NA,slope=NA,p.value=NA,p.value.negative = NA,fdr=NA,fdr.negative=NA,cell.line=cell.line)  
  }else{
    sig.df$cell.line <- cell.line
    sig.df[sig.df$fdr.negative > cut.off,'slope'] <- 0
    sig.df
  }
}



cp.across.cell.line           <- dcast(pan.can.cp.fit.df,cell.line ~ cp.name,value.var='slope') 
rownames(cp.across.cell.line) <- cp.across.cell.line$cell.line
cp.across.cell.line$cell.line <- NULL
cp.across.cell.line           <- as.matrix(cp.across.cell.line)
tmp                           <- apply(cp.across.cell.line,2, function(x) sum( x != 0 & is.na(x) == FALSE )  )
cp.across.cell.line           <- cp.across.cell.line[,order(tmp,decreasing = TRUE)]
negative.effect.cp.across.cell.line <- cp.across.cell.line


save(file='client-side/output/pan.cancer.cp.R.output/pan.cancer.cp.RData',list=c('negative.effect.cp.across.cell.line','positive.effect.cp.across.cell.line'))



# I copy and paste the compound name in negative.effect.cp.across.cell.line and positive.effect.cp.across.cell.line
# Then, I queried clue.io to get their MOA. The results were stored in client-side/meta.data/MOA.negative.hits.from.LINCS.txt and MOA.positive.hits.from.LINCS
require(data.table)
MOA.positive.hits                          <- fread(input='client-side/meta.data/MOA.positive.hits.from.LINCS.txt',header=TRUE) %>% as.data.frame
positive.effect.cp.across.cell.line.df     <- as.data.frame(positive.effect.cp.across.cell.line)
idx                                        <- match(positive.effect.cp.across.cell.line.df %>% rownames,MOA.positive.hits$Name)
positive.effect.cp.across.cell.line.df$MOA <- MOA.positive.hits$MoA[idx]
positive.effect.cp.across.cell.line.df$MOA.from.LINCS <- ifelse(is.na(positive.effect.cp.across.cell.line.df$MOA),'NO','YES')
write.csv(x=positive.effect.cp.across.cell.line.df,file='client-side/output/pan.cancer.cp.R.output/positive.effect.cp.across.cell.line.df.csv',quote=FALSE)







write.csv(x=negative.effect.cp.across.cell.line %>% t,file='client-side/output/pan.cancer.cp.R.output/negative.effect.cp.across.cell.line.csv',quote=FALSE)



# tmp <- apply(cp.across.cell.line,2, function(x) sum( is.na(x) == FALSE  ) )
# cp.across.cell.line <- cp.across.cell.line[,order(tmp,decreasing = TRUE)]
# 
# tmp <- apply(cp.across.cell.line,2, function(x) sum( is.na(x) == FALSE   ) )
# cp.across.cell.line <- cp.across.cell.line[,tmp >= 5]
# mad.vec <- apply(cp.across.cell.line,2,function(x) mad(x[is.na(x) == FALSE]))
# median.vec <- apply(cp.across.cell.line,2,function(x) median(x[is.na(x) == FALSE & x!=0  ]))





# require(plyr)
# require(dplyr)
# require(ggplot2)
# require(pheatmap)
# require(foreach)
# require(reshape2)
# 
# load('client-side/output/analyze.L1000.data.R.output/cell.line.analyzed.RData')
# 
# pan.can.oe.fit.df <- foreach(cell.line = oe.cell.line,.combine='rbind') %do% {
#   file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
#   load(file)
#   oe.fit.df$fdr.positive.effect <- p.adjust(oe.fit.df$p.value.positive.effect,method='bonferroni')
#   oe.fit.df$fdr.negative.effect <- p.adjust(oe.fit.df$p.value.negative.effect,method='bonferroni')
#   oe.fit.df$cell.line <- cell.line
#   oe.fit.df    
#   
# }
# pan.can.oe.fit.df$oe.name <- pan.can.oe.fit.df$pert_iname
# cut.off <- 0.001
# 
# sig.pan.can.positive.effect.oe.fit.df <- pan.can.oe.fit.df[pan.can.oe.fit.df$fdr.positive.effect < cut.off,c('oe.name','cell.line','p.value.positive.effect','fdr.positive.effect')]
# sig.pan.can.negative.effect.oe.fit.df <- pan.can.oe.fit.df[pan.can.oe.fit.df$fdr.negative.effect < cut.off,c('oe.name','cell.line','p.value.negative.effect','fdr.negative.effect')]
# 
# positive.cell.line.cnt.per.oe <- table(sig.pan.can.positive.effect.oe.fit.df$oe.name) %>% as.data.frame
# positive.cell.line.cnt.per.oe <- positive.cell.line.cnt.per.oe[order(positive.cell.line.cnt.per.oe$Freq,decreasing = TRUE),]
# tmp <- positive.cell.line.cnt.per.oe
# positive.cell.line.cnt.per.oe <- tmp$Freq
# names(positive.cell.line.cnt.per.oe) <- tmp$Var1
# 
# 
# 
# negative.cell.line.cnt.per.oe <- table(sig.pan.can.negative.effect.oe.fit.df$oe.name) %>% as.data.frame
# negative.cell.line.cnt.per.oe <- negative.cell.line.cnt.per.oe[order(negative.cell.line.cnt.per.oe$Freq,decreasing = TRUE),]
# tmp <- negative.cell.line.cnt.per.oe
# negative.cell.line.cnt.per.oe <- tmp$Freq
# names(negative.cell.line.cnt.per.oe) <- tmp$Var1
# 
# 
# 
# 
# 
# positive.effect.oe.across.cell.line <- dcast(sig.pan.can.positive.effect.oe.fit.df,oe.name ~ cell.line,value.var='fdr.positive.effect',fill = 1) 
# negative.effect.oe.across.cell.line <- dcast(sig.pan.can.negative.effect.oe.fit.df,oe.name ~ cell.line,value.var='fdr.negative.effect',fill = 1) 
# rownames(positive.effect.oe.across.cell.line) <- positive.effect.oe.across.cell.line$oe.name
# rownames(negative.effect.oe.across.cell.line) <- negative.effect.oe.across.cell.line$oe.name
# 
# write.csv(x=positive.effect.cp.across.cell.line ,file='client-side/output/pan.cancer.cp.R.output/positive.effect.cp.across.cell.line.csv',quote=TRUE,row.names=FALSE)
# write.csv(x=negative.effect.cp.across.cell.line ,file='client-side/output/pan.cancer.cp.R.output/negative.effect.cp.across.cell.line.csv',quote=TRUE,row.names=FALSE)
# 
# 
# 
# 
# # # 
# # # 
# # 
# # 
# # 
# # pan.can.cp.fit.df <- pan.can.cp.fit.df[complete.cases(pan.can.cp.fit.df),]
# # sig.cp.name       <- pan.can.cp.fit.df$cp.name %>% as.character %>% unique
# # 
# # pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
# #   file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
# #   load(file)
# #   sig.df <- cp.fit.df[cp.fit.df$cp.name %in% sig.cp.name,]
# #   if(sum(cp.fit.df$cp.name %in% sig.cp.name) == 0){
# #     data.frame(cp.name=NA,slope=NA,p.value=NA,p.value.negative = NA,fdr=NA,fdr.negative=NA,cell.line=cell.line)  
# #   }else{
# #     sig.df$cell.line <- cell.line
# #     sig.df[sig.df$fdr > cut.off,'slope'] <- 0
# #     sig.df
# #   }
# # }
# # 
# # 
# # 
# # cp.across.cell.line           <- dcast(pan.can.cp.fit.df,cell.line ~ cp.name,value.var='slope') 
# # rownames(cp.across.cell.line) <- cp.across.cell.line$cell.line
# # cp.across.cell.line$cell.line <- NULL
# # cp.across.cell.line           <- as.matrix(cp.across.cell.line)
# # tmp                           <- apply(cp.across.cell.line,2, function(x) sum( x != 0 & is.na(x) == FALSE )  )
# # cp.across.cell.line           <- cp.across.cell.line[,order(tmp,decreasing = TRUE)]
# # positive.effect.cp.across.cell.line <- cp.across.cell.line
# # 
# # 
# # 
# # 
# # 
# # 
# # pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
# #   file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
# #   load(file)
# #   sig.df <- cp.fit.df[cp.fit.df$fdr.negative < cut.off,]
# #   if(sum(cp.fit.df$fdr.negative < cut.off) == 0){
# #     data.frame(cp.name=NA,slope=NA,p.value=NA,p.value.negative = NA,fdr=NA,fdr.negative=NA,cell.line=cell.line)  
# #   }else{
# #     sig.df$cell.line <- cell.line
# #     sig.df
# #   }
# # }
# # pan.can.cp.fit.df <- pan.can.cp.fit.df[complete.cases(pan.can.cp.fit.df),]
# # sig.cp.name <- pan.can.cp.fit.df$cp.name %>% as.character %>% unique
# # 
# # pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
# #   file <- sprintf('client-side/output/analyze.L1000.data.R.output/%s.analyze.L1000.RData',cell.line)
# #   load(file)
# #   sig.df <- cp.fit.df[cp.fit.df$cp.name %in% sig.cp.name,]
# #   if(sum(cp.fit.df$cp.name %in% sig.cp.name) == 0){
# #     data.frame(cp.name=NA,slope=NA,p.value=NA,p.value.negative = NA,fdr=NA,fdr.negative=NA,cell.line=cell.line)  
# #   }else{
# #     sig.df$cell.line <- cell.line
# #     sig.df[sig.df$fdr.negative > cut.off,'slope'] <- 0
# #     sig.df
# #   }
# # }
# # 
# # 
# # 
# # cp.across.cell.line           <- dcast(pan.can.cp.fit.df,cell.line ~ cp.name,value.var='slope') 
# # rownames(cp.across.cell.line) <- cp.across.cell.line$cell.line
# # cp.across.cell.line$cell.line <- NULL
# # cp.across.cell.line           <- as.matrix(cp.across.cell.line)
# # tmp                           <- apply(cp.across.cell.line,2, function(x) sum( x != 0 & is.na(x) == FALSE )  )
# # cp.across.cell.line           <- cp.across.cell.line[,order(tmp,decreasing = TRUE)]
# # negative.effect.cp.across.cell.line <- cp.across.cell.line
# # 
# # 
# # save(file='client-side/output/pan.cancer.cp.R.output/pan.cancer.cp.RData',list=c('negative.effect.cp.across.cell.line','positive.effect.cp.across.cell.line'))
# # 
# # 
# # 
# # # I copy and paste the compound name in negative.effect.cp.across.cell.line and positive.effect.cp.across.cell.line
# # # Then, I queried clue.io to get their MOA. The results were stored in client-side/meta.data/MOA.negative.hits.from.LINCS.txt and MOA.positive.hits.from.LINCS
# # require(data.table)
# # MOA.positive.hits                          <- fread(input='client-side/meta.data/MOA.positive.hits.from.LINCS.txt',header=TRUE) %>% as.data.frame
# # positive.effect.cp.across.cell.line.df     <- as.data.frame(positive.effect.cp.across.cell.line)
# # idx                                        <- match(positive.effect.cp.across.cell.line.df %>% rownames,MOA.positive.hits$Name)
# # positive.effect.cp.across.cell.line.df$MOA <- MOA.positive.hits$MoA[idx]
# # positive.effect.cp.across.cell.line.df$MOA.from.LINCS <- ifelse(is.na(positive.effect.cp.across.cell.line.df$MOA),'NO','YES')
# # write.csv(x=positive.effect.cp.across.cell.line.df,file='client-side/output/pan.cancer.cp.R.output/positive.effect.cp.across.cell.line.df.csv',quote=FALSE)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # write.csv(x=negative.effect.cp.across.cell.line %>% t,file='client-side/output/pan.cancer.cp.R.output/negative.effect.cp.across.cell.line.csv',quote=FALSE)
# # 
# # 
# # 
# # # tmp <- apply(cp.across.cell.line,2, function(x) sum( is.na(x) == FALSE  ) )
# # # cp.across.cell.line <- cp.across.cell.line[,order(tmp,decreasing = TRUE)]
# # # 
# # # tmp <- apply(cp.across.cell.line,2, function(x) sum( is.na(x) == FALSE   ) )
# # # cp.across.cell.line <- cp.across.cell.line[,tmp >= 5]
# # # mad.vec <- apply(cp.across.cell.line,2,function(x) mad(x[is.na(x) == FALSE]))
# # # median.vec <- apply(cp.across.cell.line,2,function(x) median(x[is.na(x) == FALSE & x!=0  ]))
# 
