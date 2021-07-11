source('client-side/code/For.figure/ggplot.style.R')

load('client-side/output/compute.ssGSEA.score.of.signature.R.output/A375.ssGSEA.score.RData')
load('server-side/RData/lincs/A375.lincs.RData')
load('client-side/output/analyze.L1000.data.R.output/positive/A375.analyze.L1000.RData')
IFNg.ssgsea.score        <- c(cell.line.ssgsea.score.matrix['IFNg.up.gene',] - cell.line.ssgsea.score.matrix['IFNg.dn.gene',])
names(IFNg.ssgsea.score) <- colnames(cell.line.ssgsea.score.matrix)

inst.id <- cell.line.inst.info.df$inst_id[cell.line.inst.info.df$pert_iname == 'DMSO' & cell.line.inst.info.df$pert_time == 24]
plate   <- cell.line.inst.info.df$rna_plate[cell.line.inst.info.df$pert_iname == 'DMSO' & cell.line.inst.info.df$pert_time == 24]

draw.df <- data.frame(batch=plate,score=IFNg.ssgsea.score[inst.id %>% as.character],adjusted.score=adjusted.IFNg.ssgsea.score[inst.id %>% as.character])
ggplot(draw.df,aes(x=factor(plate),y=score)) + geom_boxplot() + geom_point() + ggplot.style + xlab('') + ylab('ssGSEA.score') + theme(  axis.text.x = element_blank()) + stat_summary(fun.y=median, colour="red", geom="line",group="plate",size=3)
ggplot(draw.df,aes(x=factor(plate),y=adjusted.score)) + geom_boxplot() + geom_point() + ggplot.style + xlab('') + ylab('score') + theme(  axis.text.x = element_blank())

inst.id <- cell.line.inst.info.df$inst_id[ cell.line.inst.info.df$pert_time == 24 & cell.line.inst.info.df$pert_type == 'trt_cp']
plate   <- cell.line.inst.info.df$rna_plate[ cell.line.inst.info.df$pert_time == 24 & cell.line.inst.info.df$pert_type == 'trt_cp']
draw.df <- data.frame(batch=plate,score=IFNg.ssgsea.score[inst.id %>% as.character],adjusted.score=adjusted.IFNg.ssgsea.score[inst.id %>% as.character])
ggplot(draw.df,aes(x=factor(plate),y=score)) + geom_boxplot() + ggplot.style + xlab('') + ylab('ssGSEA.score') + theme(  axis.text.x = element_blank()) + stat_summary(fun.y=median, colour="red", geom="line",group="plate",size=3)
ggplot(draw.df,aes(x=factor(plate),y=adjusted.score)) + geom_boxplot() +  ggplot.style + xlab('') + ylab('ssGSEA.score') + theme(  axis.text.x = element_blank())+ stat_summary(fun.y=median, colour="red", geom="line",group="plate",size=3)
