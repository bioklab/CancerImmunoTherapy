# This is the main file implementing the workflow!
require(GSVA)
require(data.table)
require(plyr)
require(dplyr)
require(foreach)
require(ggplot2)
require(quantreg)

combine.p.value.fisher <- function(df) {
    fisher.statistics     <- -2 * sum(log(df$p.value))
    p.value               <- pchisq(q=fisher.statistics,df = 2 * nrow(df),lower.tail = FALSE)
    rs                    <- data.frame(p.value = p.value,pert.number=nrow(df))
    
    if(nrow(df) == 1){
        rs$p.value        <- df$p.value 
        rs$effect.size    <- df$adjusted.score / df$pert_dose
    }
    else{
        lm.rs          <- lm(data = df,formula =  adjusted.score ~ pert_dose + 0)
        rs$effect.size <- lm.rs$coefficients[1]
    }
    rs
}


df            <- read.csv(file='server-side/RData/lincs/LINCS.data.statistics.csv',header=TRUE)
cp.cell.line  <- df$cell_id[df$trt_cp > 6000] %>% as.character
oe.cell.line  <- df$cell_id[df$trt_oe > 0]    %>% as.character
#cell.line.vec <- c(cp.cell.line,oe.cell.line) %>% unique

for(cell.line in oe.cell.line){
  
  lincs.data.file   <- sprintf('server-side/RData/lincs/%s.lincs.RData',cell.line)
  ssgsea.score.file <- sprintf('client-side/output/compute.ssGSEA.score.of.signature.R.output/%s.ssGSEA.score.RData',cell.line)
  rdata.file        <- sprintf('client-side/output/analyze.L1000.data.R.output/positive/%s.analyze.L1000.RData',cell.line)
  if(file.exists(rdata.file)){
    next
  }
  
  load(lincs.data.file)
  load(ssgsea.score.file)
  
  
  rownames(cell.line.inst.info.df) <- cell.line.inst.info.df$inst_id
  
  #IFNg.ssgsea.score        <- c(cell.line.ssgsea.score.matrix['IFN.gamma',]) # this is MSigDB signature ssGSEA score
  IFNg.ssgsea.score        <- c(cell.line.ssgsea.score.matrix['IFNg.dn.gene',]- cell.line.ssgsea.score.matrix['IFNg.up.gene',])
  
  names(IFNg.ssgsea.score) <- colnames(cell.line.ssgsea.score.matrix)
  
  
  
  ################# estimation of plate effect ###############
  rna.plate.vec        <- cell.line.inst.info.df$rna_plate %>% as.character %>% unique
  f                    <- grepl(x=cell.line.inst.info.df$pert_type,pattern = 'ctl')
  ctrl.pert.vec        <- cell.line.inst.info.df$pert_iname[f] %>% unique %>% as.character
  plate.effect.matrix  <- foreach(rna.plate = rna.plate.vec,.combine='rbind') %do% {
    flag    <- cell.line.inst.info.df$rna_plate == rna.plate 
    inst.id <- cell.line.inst.info.df$inst_id[flag] %>% as.character
    e1      <- median(IFNg.ssgsea.score[inst.id]) # median across all perturbagens
    
    v <- foreach(ctrl.pert = ctrl.pert.vec,.combine='c') %do% {
      flag    <- cell.line.inst.info.df$rna_plate == rna.plate & cell.line.inst.info.df$pert_iname == ctrl.pert
      if(sum(flag) == 0){    
        e2 <- NA
      }else{
        inst.id <- cell.line.inst.info.df$inst_id[flag] %>% as.character
        e2      <- median(IFNg.ssgsea.score[inst.id]) # median across control  perturbagens
      }
      e2
    }
    c(e1,v)    
  }
  rownames(plate.effect.matrix) <- rna.plate.vec
  colnames(plate.effect.matrix) <- c('median.of.all.pert',ctrl.pert.vec)
  
  
  trt.oe.plate  <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_oe','rna_plate']  %>% unique %>% as.character
  trt.sh.plate  <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_sh','rna_plate']  %>% unique %>% as.character
  trt.cp.plate  <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_cp','rna_plate']  %>% unique %>% as.character
  trt.lig.plate <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_lig','rna_plate'] %>% unique %>% as.character
  
  
  plate.effect.vec <- plate.effect.matrix[,'median.of.all.pert']
  plate.effect.vec[trt.cp.plate] <- plate.effect.matrix[trt.cp.plate,'DMSO']
  if(length(trt.oe.plate) > 0){
      plate.effect.vec[trt.oe.plate] <- plate.effect.matrix[trt.oe.plate,'UnTrt']
  }
  
  
  idx              <- cell.line.inst.info.df$rna_plate[match(names(IFNg.ssgsea.score),cell.line.inst.info.df$inst_id)] %>% as.character
  adjusted.IFNg.ssgsea.score <- IFNg.ssgsea.score - plate.effect.vec[idx]
  
  
  ################################### ranking ligand  pert##############
  lig.df       <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_lig' & cell.line.inst.info.df$pert_time == 4,]
  lig.df$score <- adjusted.IFNg.ssgsea.score[rownames(lig.df)]
  tmp          <- ddply(lig.df,.(pert_iname),function(x) data.frame(score=median(x$score),mad=mad(x$score)))
  lig.fit.df   <- tmp
  lig.fit.df   <- lig.fit.df[order(lig.fit.df$score,decreasing = TRUE),]
    
    
  
  
  
  
  
  
  ################################### ranking overepxression pert ##############
  
  oe.plate    <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_oe' & cell.line.inst.info.df$pert_time == 96,'rna_plate'] %>% unique %>% as.character
  if(cell.line == 'HEK293T'){
    oe.plate    <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_oe' & cell.line.inst.info.df$pert_time == 48,'rna_plate'] %>% unique %>% as.character
  }
  if(length(oe.plate) > 0){
    oe.df       <- cell.line.inst.info.df[cell.line.inst.info.df$rna_plate %in% oe.plate,]
    oe.df$adjusted.score <- adjusted.IFNg.ssgsea.score[rownames(oe.df)]
    
    pert.num.df <- table(oe.df$pert_iname) %>% as.data.frame
    pert.num.df <- pert.num.df[order(pert.num.df$Freq,decreasing = TRUE),]
    
    #sigma.estimate <- 1.4826 * mad(adjusted.IFNg.ssgsea.score[oe.df$inst_id[oe.df$pert_iname == pert.num.df$Var1[1]] %>% as.character]) 
    sigma.estimate <- 1.4826 * mad(adjusted.IFNg.ssgsea.score[oe.df$inst_id[oe.df$pert_iname == 'UnTrt'] %>% as.character]) 
    
    oe.df$p.value <- pnorm(q = oe.df$adjusted.score,mean = 0,sd=sigma.estimate,lower.tail = FALSE)
    #oe.df$p.value.negative.effect <- pnorm(q = oe.df$score,mean = 0,sd=sigma.estimate,lower.tail = TRUE)
    
    tmp           <- ddply(oe.df,.(pert_iname),
                           function(x) {
                             test.statistics               <- -2 * sum(log(x$p.value))
                             p.value       <- pchisq(q=test.statistics,df = 2 * nrow(x),lower.tail = FALSE) 
                             data.frame(p.value=p.value,effect.size=median(x$adjusted.score))
                           }
    )
    tmp                 <- tmp[order(tmp$p.value),]
    oe.fit.df           <- tmp
    rownames(oe.fit.df) <- oe.fit.df$pert_iname
  }else{
    oe.fit.df <- NULL  
  }
  
  
  ###########################
  
  save(file=rdata.file,list=c('oe.fit.df','lig.fit.df','adjusted.IFNg.ssgsea.score','plate.effect.matrix','trt.oe.plate','trt.cp.plate','trt.sh.plate','trt.lig.plate'))
  
}

oe.fit.df$fdr <- p.adjust(oe.fit.df$p.value,method='fdr')
sig.df        <- oe.fit.df[oe.fit.df$fdr < 0.05,]
sig.df        <- sig.df[order(sig.df$effect.size,decreasing = TRUE),]

####################################################### shRNA pertrubagen ###################################
for(cell.line in cell.line.vec){
  
  lincs.data.file   <- sprintf('server-side/RData/lincs/%s.lincs.RData',cell.line)
  ssgsea.score.file <- sprintf('client-side/output/compute.ssGSEA.score.of.signature.R.output/%s.ssGSEA.score.RData',cell.line)
  rdata.file        <- sprintf('client-side/output/analyze.L1000.data.R.output/negative/%s.analyze.L1000.RData',cell.line)
  if(file.exists(rdata.file)){
    next
  }
  
  load(lincs.data.file)
  load(ssgsea.score.file)
  
  
  rownames(cell.line.inst.info.df) <- cell.line.inst.info.df$inst_id
  
  #IFNg.ssgsea.score        <- c(cell.line.ssgsea.score.matrix['IFN.gamma',]) # this is MSigDB signature ssGSEA score
  IFNg.ssgsea.score        <- c( cell.line.ssgsea.score.matrix['IFNg.up.gene',] - cell.line.ssgsea.score.matrix['IFNg.dn.gene',])
  IFNg.ssgsea.score        <- c( cell.line.ssgsea.score.matrix['IFNg.up.gene',] )
  
  names(IFNg.ssgsea.score) <- colnames(cell.line.ssgsea.score.matrix)
  
  
  
  ################# estimation of plate effect ###############
  rna.plate.vec        <- cell.line.inst.info.df$rna_plate %>% as.character %>% unique
  f                    <- grepl(x=cell.line.inst.info.df$pert_type,pattern = 'ctl')
  ctrl.pert.vec        <- cell.line.inst.info.df$pert_iname[f] %>% unique %>% as.character
  plate.effect.matrix  <- foreach(rna.plate = rna.plate.vec,.combine='rbind') %do% {
    flag    <- cell.line.inst.info.df$rna_plate == rna.plate 
    inst.id <- cell.line.inst.info.df$inst_id[flag] %>% as.character
    e1      <- median(IFNg.ssgsea.score[inst.id]) # median across all perturbagens
    
    v <- foreach(ctrl.pert = ctrl.pert.vec,.combine='c') %do% {
      flag    <- cell.line.inst.info.df$rna_plate == rna.plate & cell.line.inst.info.df$pert_iname == ctrl.pert
      if(sum(flag) == 0){    
        e2 <- NA
      }else{
        inst.id <- cell.line.inst.info.df$inst_id[flag] %>% as.character
        e2      <- median(IFNg.ssgsea.score[inst.id]) # median across control  perturbagens
      }
      e2
    }
    c(e1,v)    
  }
  rownames(plate.effect.matrix) <- rna.plate.vec
  colnames(plate.effect.matrix) <- c('median.of.all.pert',ctrl.pert.vec)
  
  
  trt.oe.plate  <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_oe','rna_plate']  %>% unique %>% as.character
  trt.sh.plate  <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_sh','rna_plate']  %>% unique %>% as.character
  trt.cp.plate  <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_cp','rna_plate']  %>% unique %>% as.character
  trt.lig.plate <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type =='trt_lig','rna_plate'] %>% unique %>% as.character
  
  
  plate.effect.vec <- plate.effect.matrix[,'median.of.all.pert']
  plate.effect.vec[trt.cp.plate] <- plate.effect.matrix[trt.cp.plate,'DMSO']
  if(length(trt.oe.plate) > 0){
    plate.effect.vec[trt.oe.plate] <- plate.effect.matrix[trt.oe.plate,'UnTrt']
  }
  if(length(trt.sh.plate) > 0){
    plate.effect.vec[trt.sh.plate] <- plate.effect.matrix[trt.sh.plate,'GFP']
  }
  
  
  idx              <- cell.line.inst.info.df$rna_plate[match(names(IFNg.ssgsea.score),cell.line.inst.info.df$inst_id)] %>% as.character
  adjusted.IFNg.ssgsea.score <- IFNg.ssgsea.score - plate.effect.vec[idx]
  
  
  
  
  
  ################################### ranking ligand  pert##############
  lig.df       <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_lig' & cell.line.inst.info.df$pert_time == 4,]
  lig.df$score <- adjusted.IFNg.ssgsea.score[rownames(lig.df)]
  tmp          <- ddply(lig.df,.(pert_iname),function(x) data.frame(score=median(x$score),mad=mad(x$score)))
  lig.fit.df   <- tmp
  lig.fit.df   <- lig.fit.df[order(lig.fit.df$score,decreasing = TRUE),]
  
  
  
  
  ################################### ranking overepxression pert ##############
  
  sh.plate    <- cell.line.inst.info.df[cell.line.inst.info.df$pert_type == 'trt_sh','rna_plate'] %>% unique %>% as.character
  if(length(sh.plate) > 0){
    sh.df       <- cell.line.inst.info.df[cell.line.inst.info.df$rna_plate %in% sh.plate,]
    sh.df$adjusted.score <- adjusted.IFNg.ssgsea.score[rownames(sh.df)]
    
    pert.num.df <- table(sh.df$pert_iname) %>% as.data.frame
    pert.num.df <- pert.num.df[order(pert.num.df$Freq,decreasing = TRUE),]
    
    #sigma.estimate <- 1.4826 * mad(adjusted.IFNg.ssgsea.score[sh.df$inst_id[sh.df$pert_iname == pert.num.df$Var1[1]] %>% as.character]) 
    sigma.estimate <- 1.4826 * mad(adjusted.IFNg.ssgsea.score[sh.df$inst_id[sh.df$pert_iname == 'UnTrt'] %>% as.character]) 
    
    sh.df$p.value <- pnorm(q = sh.df$adjusted.score,mean = 0,sd=sigma.estimate,lower.tail = FALSE)
    #sh.df$p.value.negative.effect <- pnorm(q = sh.df$score,mean = 0,sd=sigma.estimate,lower.tail = TRUE)
    
    tmp           <- ddply(sh.df,.(pert_iname),
                           function(x) {
                             test.statistics               <- -2 * sum(log(x$p.value))
                             p.value       <- pchisq(q=test.statistics,df = 2 * nrow(x),lower.tail = FALSE) 
                             data.frame(p.value=p.value,effect.size=median(x$adjusted.score))
                           }
    )
    tmp                 <- tmp[order(tmp$p.value),]
    sh.fit.df           <- tmp
    rownames(sh.fit.df) <- sh.fit.df$pert_iname
  }else{
    sh.fit.df <- NULL  
  }
  
  
  ###########################
  
  
}



sh.fit.df$fdr <- p.adjust(sh.fit.df$p.value,method='fdr')
sig.df        <- sh.fit.df[sh.fit.df$fdr < 0.05,]
sig.df        <- sig.df[order(sig.df$effect.size,decreasing = TRUE),]
View(sig.df)













############################################### garbage code #####################################

# estimate.sigma <- function(cp.df,pert.time){
#     tmp <- cp.df[cp.df$pert_time == pert.time ,]     
#     rs <- ddply(tmp,.(pert_iname,pert_time,pert_dose),nrow)  
#     rs <- rs[order(rs$V1,decreasing = TRUE),]
#     pert.iname <- rs$pert_iname[1]
#     pert.dose <- rs$pert_dose[1]
#     inst.id <- cp.df$inst_id[cp.df$pert_time == pert.time & cp.df$pert_iname == pert.iname & cp.df$pert_dose == pert.dose]
#     print(sprintf('%s %s',pert.time,pert.iname))
#     1.4826 * mad(adjusted.IFNg.ssgsea.score[inst.id])
# }



# A375.effect.size.vec <- effect.size.vec
# A375.p.value.vec <- p.value.vec
# 
# A549.effect.size.vec <- effect.size.vec
# A549.p.value.vec <- p.value.vec
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# df <- data.frame(fdr=-1 * (fdr.vec[common.cp] %>% log10),effect.size=effect.size.vec[common.cp])
# df$cp.name <- rownames(df)
# require(plotly)
# 
# 
# ggplot(df,aes(x=effect.size,y=fdr,label=cp.name)) + geom_text() 
# 
# pert.str                   <- cell.line.inst.info.df$pert_str[match(names(adjusted.IFNg.ssgsea.score),cell.line.inst.info.df$inst_id)]
# flag1 <- pert.str == 'DMSO_ctl_vehicle_-666_24'
# flag2 <- pert.str == 'DMSO_ctl_vehicle_-666_6'
# flag3 <- pert.str == 'vorinostat_trt_cp_0.1_24'
# flag4 <- pert.str == 'vorinostat_trt_cp_10_24'
# 
# var.test(adjusted.IFNg.ssgsea.score[flag3],adjusted.IFNg.ssgsea.score[flag4],alternative = 'less')




# mad.df <- ddply(cp.df,.(rna_plate),function(x) mad(x$score))
# mad.df <- mad.df[order(mad.df$V1),]
# cp.df <- merge(cp.df,mad.df,by='rna_plate')
# View(cp.df[order(cp.df$score),])
# 
# rrs <- foreach(r= cell.line.inst.info.df$rna_plate %>% unique %>% as.character,.combine='rbind') %do% {
# flag <- cell.line.inst.info.df$rna_plate == r
# cp.df <- cell.line.inst.info.df[flag,]
# cp.df$score  <- adjusted.IFNg.ssgsea.score[rownames(cp.df) %>% as.character]
# df <- ddply(cp.df,.(pert_str), function(x) data.frame(m=median(x$score),v=sd(x$score),nu=nrow(x)) )
# #plot(x=df$m,y=df$v)
# df[order(df$v),] 
# }







# target.pert     <- pert.num.df$pert_str[pert.num.df$pert.num >= 3] %>% as.character
# pert.str        <- cell.line.inst.info.df$pert_str[match(names(adjusted.IFNg.ssgsea.score),cell.line.inst.info.df$inst_id)]
# flag            <- pert.str %in% target.pert
# combined.score  <- tapply(adjusted.IFNg.ssgsea.score[flag],pert.str[flag],median)
# combined.score  <- c(combined.score)
# combined.score  <- sort(combined.score,decreasing = TRUE)
# 
# 
# 
# y  <- tapply(IFNg.ssgsea.score[flag],pert.str[flag],median)
# y  <- c(y)
# y  <- sort(y,decreasing = TRUE)


# tmp <- ddply(cell.line.inst.info.df,.(rna_plate),nrow)
# View(tmp)
# 
# 
# common.cp <- intersect(names(effect.size.vec),names(p.value.vec))
# plot(x=effect.size.vec[common.cp],y=-1 * (fdr.vec[common.cp] %>% log10))
# abline(h=3)
# df <- foreach(str = tmp$pert_str[tmp$pert.num >= 3],.combine='rbind') %do% {
#     inst.id <- cell.line.inst.info.df[cell.line.inst.info.df$pert_str == str,'inst_id'] %>% as.character
#     data.frame(median=median(adjusted.IFNg.ssgsea.score[inst.id]),mad=mad(adjusted.IFNg.ssgsea.score[inst.id]))
# }
# plot(x=df$median,y=df$mad)
# rownames(df) <- tmp$pert_str[tmp$pert.num >= 3]
# df <- df[order(df$median),]
# 
# cp.df <- df[grepl(x=rownames(df),pattern = 'trt_cp'),]
# cp.df <- cp.df[order(cp.df$mad,decreasing = TRUE),]
# 


# control.perturbation <- c('DMSO','UnTrt','lacZ','LUCIFERASE','EMPTY_VECTOR','GFP','pgw','RFP','HcRed','eGFP')
# control.perturbation <- c('DMSO')
#     inst.id1 <- cell.line.inst.info.df$inst_id[flag1]
#     flag2    <- (cell.line.inst.info.df$rna_plate == rna.plate) & (cell.line.inst.info.df$pert_iname %in% control.perturbation)
#     inst.id2 <- cell.line.inst.info.df$inst_id[flag2]
#all.effect  <- ifelse( length(inst.id1) == 0, -1024,median(IFNg.ssgsea.score[inst.id1]) )
#ctrl.effect <- ifelse( length(inst.id2) == 0, -1024,median(IFNg.ssgsea.score[inst.id2]) )
#c(all.effect,ctrl.effect)
# common.gene <- intersect(L1000.gene.entrez.id,genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']])
# median.expr <- apply(cell.line.inst.matrix[common.gene,],2,median)
# 
# idx <- sample(x = 1:10000,size=5000)
# plot(x=median.expr[idx],y=IFNg.ssgsea.score[idx])
# tt <- cell.line.inst.info.df[,c('inst_id','rna_plate')]
# tt$score <- IFNg.ssgsea.score[tt[,'inst_id'] %>% as.character]
# 
# plate.score <- ddply(tt,.(rna_plate),function(x) data.frame(median=mean(x$score),mad=mad(x$score)))
# 
# ttex <- merge(tt,plate.score)
# 
# ttex$new.score <- ttex$score - ttex$median
# IFNg.ssgsea.new.score <- ttex$new.score
# names(IFNg.ssgsea.new.score) <- (ttex$inst_id)
# 
# 
# df <- foreach(str = tmp$pert_str[tmp$pert.num >= 3],.combine='rbind') %do% {
#   inst.id <-   cell.line.inst.info.df[cell.line.inst.info.df$pert_str == str,'inst_id'] %>% as.character
#   data.frame(median=median(IFNg.ssgsea.new.score[inst.id]),mad=mad(IFNg.ssgsea.new.score[inst.id]))
# }
# plot(x=df$median,y=df$mad)
# rownames(df) <- tmp$pert_str[tmp$pert.num >= 3]
# df <- df[order(df$median),]
# 
# cp.df <- df[grepl(x=rownames(df),pattern = 'trt_cp'),]
# cp.df <- cp.df[order(cp.df$median,decreasing = TRUE),]
# 
# dose <- sapply(strsplit(x = rownames(cp.df),split='_'),function(x) x[[4]])
# time <- sapply(strsplit(x = rownames(cp.df),split='_'),function(x) x[[5]])
# drug <- sapply(strsplit(x = rownames(cp.df),split='_'),function(x) x[[1]])
# 
# flag <- drug == 'curcumin' & time == '24'
# 
# plot(x=dose[flag],y=cp.df$median[flag])
# 
# 
# sh.df <- df[grepl(x=rownames(df),pattern = 'trt_sh'),]
# sh.df <- sh.df[order(sh.df$median,decreasing = TRUE),]

# r <- apply(cell.line.inst.matrix,2,rank)
# s <- apply(r,1,mad)
# d <- apply(r,1,median)
# 
# 
# cmpd.up.gene.score <- c( t(cmpd.expr[sig.up.gene,]) %*% (-1 * log10(rna.seq.signature.matrix[sig.up.gene,'fdr'])) )
# cmpd.up.gene.score <- c( t(cmpd.expr[sig.up.gene,]) %*% ((rna.seq.signature.matrix[sig.up.gene,'adjusted.cor'])) )
# 
# names(cmpd.up.gene.score) <- colnames(cmpd.expr)
# cmpd.up.gene.score         <- sort(cmpd.up.gene.score,decreasing = TRUE)
# signature_meta[names(cmpd.up.gene.score) %>% head(300),] %>% View
# cmpd.name <- signature_meta[names(cmpd.up.gene.score) ,'pert_iname']

# plate.vec <- cell.line.inst.info.df$rna_plate %>% unique  %>% as.character
# f.test.p.value.df <- foreach(rna.plate = plate.vec[150:160],.combine='rbind') %do% {
#     flag1      <- cell.line.inst.info.df$rna_plate == rna.plate &cell.line.inst.info.df$pert_iname == 'DMSO'
#     if(sum(flag1) ==0){return(data.frame(p.value=NA,cp.name=NA,pert.time=NA))}
#     inst.id1   <- cell.line.inst.info.df$inst_id[flag1] %>% as.character
#     dmso.score <- IFNg.ssgsea.score[inst.id1]
#     cp.name    <- cell.line.inst.info.df[cell.line.inst.info.df$rna_plate == rna.plate,'pert_iname'] %>% unique %>% as.character
#     v <- foreach(cp = cp.name,.combine='c') %do% {
#         flag2      <- cell.line.inst.info.df$rna_plate == rna.plate & cell.line.inst.info.df$pert_iname == cp
#         inst.id2   <- cell.line.inst.info.df$inst_id[flag2] %>% as.character
#         if(length(inst.id2) == 1) {return(NA)}
#         var.test(dmso.score,IFNg.ssgsea.score[inst.id2],alternative = 'less')$p.value
#     }
#     data.frame(p.value=v,cp.name=cp.name,pert.time=cell.line.inst.info.df$pert_time[flag1] %>% unique)
# }
# f.test.p.value.df <- f.test.p.value.df[order(f.test.p.value.df$p.value),]
# View(f.test.p.value.df)
# 
# 
# 
# 
# gel.inst.id <- cell.line.inst.info.df[cell.line.inst.info.df$pert_str == 'DMSO_ctl_vehicle_-666_24','inst_id']  %>% as.character
# df <- data.frame(score=IFNg.ssgsea.score[gel.inst.id],plate=cell.line.inst.info.df[gel.inst.id,'rna_plate'])
# #ggplot(df) + geom_boxplot(aes(x=plate,y=score))
# medin.df <- ddply(df,.(plate),function(x) median(x$score))
# 
# 
# 
# mad.df <- ddply(df,.(plate),function(x) mad(x$score))
# rownames(mad.df) <- mad.df$plate
# 
# df <- data.frame(score=IFNg.ssgsea.score,plate=cell.line.inst.info.df[names(IFNg.ssgsea.score),'rna_plate'])
# #ggplot(df) + geom_boxplot(aes(x=plate,y=score))
# mad.plate.df <- ddply(df,.(plate),function(x) mad(x$score))
# rownames(mad.plate.df) <- mad.plate.df$plate
# 
# plot(x=mad.df[,'V1'],y=mad.plate.df[rownames(mad.df),'V1'],xlab='DMSO var',ylab='plate var',xlim=c(0,30),ylim=c(0,30))
# lines(c(0,30),c(0,30))
# 
# hehe <- mad.plate.df[rownames(mad.df),'V1'] * mad.plate.df[rownames(mad.df),'V1'] - mad.df[,'V1'] * mad.df[,'V1']
# 
# time.vec <- sapply( strsplit(x=mad.plate.df$plate %>% as.character,split = '_'), function(x) x[[3]])
# mad.plate.df$time <- time.vec
# ggplot(mad.plate.df,aes(x=time,y=V1)) + geom_boxplot()
# names(p.value.vec) <- cp.name
# 
# p.value.vec <- p.value.vec[p.value.vec >0]
# hist(p.value.vec)
# View(sort(p.value.vec))
# fdr.vec <- p.adjust(p.value.vec,method='fdr')
# 
# flag <- cell.line.inst.info.df$pert_iname == 'genz-644282' & cell.line.inst.info.df$pert_time == 24
# cp.df <- cell.line.inst.info.df[flag,]
# cp.df$score  <- adjusted.IFNg.ssgsea.score[rownames(cp.df) %>% as.character]
# ggplot(rbind(cp.df[,c('score','pert_dose')],zero.df),aes(x=pert_dose,y=score)) + geom_point()  + stat_smooth(method='lm')
# View(cp.df[order(cp.df$score),])
# fit.df <- fit.df[order(fit.df$slope,decreasing = TRUE),]
# rownames(fit.df) <- fit.df$cp.name
# 
# 
# r1 <- which(A375.cp.lincs.query.df$Name %in% HDACi) / nrow(A375.cp.lincs.query.df)
# r2 <- which(fit.df$cp.name %in% HDACi) / nrow(fit.df)
# 
# common.cp<- intersect(A375.cp.lincs.query.df$Name,fit.df$cp.name)
# 
# A375.cp.lincs.query.df <- A375.cp.lincs.query.df[A375.cp.lincs.query.df$Name %in% common.cp,]
# haha <- fit.df[fit.df$cp.name %in% common.cp,]
# 
# 
# hc <- intersect(common.cp,HDACi)
# hc <- intersect(common.cp,CDKi)
# 
# A375.cp.lincs.query.df$rank <- 1:nrow(A375.cp.lincs.query.df)
# haha$rank <- 1:nrow(haha)
# 
# r1 <- A375.cp.lincs.query.df[hc,'rank']
# r2 <- haha[hc,'rank']
# 
# 
# effect.size.vec <- foreach(cp = cp.name %>% as.character,.combine='c') %do% {
#   cp.inst.id <- rownames(cp.df)[cp.df$pert_iname == cp]
#   #if(length(cp.inst.id) <= 10){return(-1)}
#   median(adjusted.IFNg.ssgsea.score[cp.inst.id]) - median(adjusted.IFNg.ssgsea.score[dmso.inst.id]) 
# }
# names(effect.size.vec) <- cp.name
#df                <- fread(input = 'client-side/meta.data/L1000.land.mark.gene.tsv.txt') %>% as.data.frame
#land.mark.gene.id <- df[,1] %>% as.character
#load('server-side/RData/lincs/A375.lincs.RData')
#load('client-side/output/compute.ssGSEA.score.of.signature.R.output/A375.ssGSEA.score.RData')
#load('client-side//output//prepare.meta.data.R.output//prepare.meta.data.RData')
# cell.line.inst.info.df$pert_str  <- paste(cell.line.inst.info.df$pert_iname, cell.line.inst.info.df$pert_type,cell.line.inst.info.df$pert_dose,cell.line.inst.info.df$pert_time,sep='_')
# tmp                              <- ddply(cell.line.inst.info.df,.(pert_str),nrow)
# tmp                              <- tmp[order(tmp$V1,decreasing = TRUE),]
# colnames(tmp)[2]                 <- 'pert.num'
# pert.num.df                      <- tmp
