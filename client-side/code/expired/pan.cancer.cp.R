require(plyr)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(foreach)
require(reshape2)

load('client-side/output/analyze.L1000.data.R.output/cell.line.analyzed.RData')
#positive cp
pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
    file         <- sprintf('client-side/output/analyze.L1000.data.R.output/positive/%s.analyze.L1000.RData',cell.line)
    load(file)
    cp.fit.df$fdr <- p.adjust(cp.fit.df$p.value,method='bonferroni')
    cp.fit.df$cell.line <- cell.line
    cp.fit.df    
    
}
pan.can.cp.fit.df$cp.name <- pan.can.cp.fit.df$pert_iname
cut.off                   <- 0.001

pan.can.cp.fit.wide.df <- dcast(pan.can.cp.fit.df,cp.name ~ cell.line,value.var='fdr')
sig.cp                 <- pan.can.cp.fit.df$cp.name[pan.can.cp.fit.df$fdr < cut.off] %>% unique %>% as.character # significant in at least one cell line
sig.cp.across.cell.line           <- pan.can.cp.fit.wide.df[pan.can.cp.fit.wide.df$cp.name %in% sig.cp,]
rownames(sig.cp.across.cell.line) <- sig.cp.across.cell.line$cp.name
sig.cp.across.cell.line$cp.name   <- NULL
sig.cp.across.cell.line           <- as.matrix(sig.cp.across.cell.line)
sig.cell.line.cnt.per.cp          <- apply(sig.cp.across.cell.line,1,function(x){
    y <- x[is.na(x) == FALSE]
    sum(y < cut.off)
})
sig.cell.line.cnt.per.cp          <- sort(sig.cell.line.cnt.per.cp,decreasing = TRUE)
positive.sig.cell.line.cnt.per.cp <- sig.cell.line.cnt.per.cp
positive.sig.cp.acorss.cell.line  <- sig.cp.across.cell.line


### negative cp
pan.can.cp.fit.df <- foreach(cell.line = cp.cell.line,.combine='rbind') %do% {
  file         <- sprintf('client-side/output/analyze.L1000.data.R.output/negative/%s.analyze.L1000.RData',cell.line)
  load(file)
  cp.fit.df$fdr <- p.adjust(cp.fit.df$p.value,method='bonferroni')
  cp.fit.df$cell.line <- cell.line
  cp.fit.df    
  
}
pan.can.cp.fit.df$cp.name <- pan.can.cp.fit.df$pert_iname
cut.off                   <- 0.001

pan.can.cp.fit.wide.df <- dcast(pan.can.cp.fit.df,cp.name ~ cell.line,value.var='fdr')
sig.cp                 <- pan.can.cp.fit.df$cp.name[pan.can.cp.fit.df$fdr < cut.off] %>% unique %>% as.character # significant in at least one cell line
sig.cp.across.cell.line           <- pan.can.cp.fit.wide.df[pan.can.cp.fit.wide.df$cp.name %in% sig.cp,]
rownames(sig.cp.across.cell.line) <- sig.cp.across.cell.line$cp.name
sig.cp.across.cell.line$cp.name   <- NULL
sig.cp.across.cell.line           <- as.matrix(sig.cp.across.cell.line)
sig.cell.line.cnt.per.cp          <- apply(sig.cp.across.cell.line,1,function(x){
  y <- x[is.na(x) == FALSE]
  sum(y < cut.off)
})
sig.cell.line.cnt.per.cp          <- sort(sig.cell.line.cnt.per.cp,decreasing = TRUE)
negative.sig.cell.line.cnt.per.cp <- sig.cell.line.cnt.per.cp
negative.sig.cp.acorss.cell.line  <- sig.cp.across.cell.line

save(file='client-side/output/pan.cancer.cp.R.output/pan.cancer.cp.RData',list=c('negative.sig.cell.line.cnt.per.cp','positive.sig.cell.line.cnt.per.cp'))






#write.csv(x=positive.effect.cp.across.cell.line ,file='client-side/output/pan.cancer.cp.R.output/positive.effect.cp.across.cell.line.csv',quote=TRUE,row.names=FALSE)
#write.csv(x=negative.effect.cp.across.cell.line ,file='client-side/output/pan.cancer.cp.R.output/negative.effect.cp.across.cell.line.csv',quote=TRUE,row.names=FALSE)

positive.cp <- positive.sig.cell.line.cnt.per.cp %>% names
negative.cp <- negative.sig.cell.line.cnt.per.cp %>% names
intersect(positive.cp,negative.cp)


rownames(positive.effect.cp.across.cell.line) <- positive.effect.cp.across.cell.line$cp.name 
rownames(negative.effect.cp.across.cell.line) <- negative.effect.cp.across.cell.line$cp.name 

positive.effect.cp.across.cell.line$cp.name <- NULL
negative.effect.cp.across.cell.line$cp.name <- NULL

positive.effect.cp.across.cell.line <- as.matrix(positive.effect.cp.across.cell.line)
negative.effect.cp.across.cell.line <- as.matrix(negative.effect.cp.across.cell.line)


negative.cell.line.cnt.per.cp <- apply(negative.effect.cp.across.cell.line,1,function(x){
  y <- x[is.na(x) == FALSE]
  sum(y < cut.off)
})
negative.cell.line.cnt.per.cp <- sort(negative.cell.line.cnt.per.cp,decreasing = TRUE)



positive.cell.line.cnt.per.cp <- apply(positive.effect.cp.across.cell.line,1,function(x){
  y <- x[is.na(x) == FALSE]
  sum(y < cut.off)
})
positive.cell.line.cnt.per.cp <- sort(positive.cell.line.cnt.per.cp,decreasing = TRUE)


HDACi <- read.table('client-side/external.data/from.LINCS//LINCS.HDACi.txt',header=FALSE)$V1 %>% as.character %>% unique
CDKi  <- read.table('client-side/external.data/from.LINCS/LINCS.CDKi.txt',header=FALSE)$V1 %>% as.character  %>% unique
RAFi  <- read.table('client-side/external.data/from.LINCS/LINCS.RAFi.txt',header=FALSE)$V1 %>% as.character  %>% unique
HSPi  <- read.table('client-side/external.data/from.LINCS/LINCS.HSPi.txt',header=FALSE)$V1 %>% as.character  %>% unique
DNMTi  <- read.table('client-side/external.data/from.LINCS/LINCS.DNMTi.txt',header=FALSE)$V1 %>% as.character  %>% unique


positive.cell.line.cnt[CDKi]
positive.cell.line.cnt[HDACi]
positive.cell.line.cnt[RAFi]




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

