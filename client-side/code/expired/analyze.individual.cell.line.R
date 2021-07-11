
require(data.table)

HDACi <- read.table('client-side/external.data/from.LINCS//LINCS.HDACi.txt',header=FALSE)$V1 %>% as.character %>% unique
CDKi  <- read.table('client-side/external.data/from.LINCS/LINCS.CDKi.txt',header=FALSE)$V1 %>% as.character  %>% unique
RAFi  <- read.table('client-side/external.data/from.LINCS/LINCS.RAFi.txt',header=FALSE)$V1 %>% as.character  %>% unique
HSPi  <- read.table('client-side/external.data/from.LINCS/LINCS.HSPi.txt',header=FALSE)$V1 %>% as.character  %>% unique
DNMTi  <- read.table('client-side/external.data/from.LINCS/LINCS.DNMTi.txt',header=FALSE)$V1 %>% as.character  %>% unique
cut.off <- 0.001
MOA.df <- fread(input='client-side/external.data/from.LINCS/TOUCHSTONE.compound.MOA.txt',header=TRUE) %>% as.data.frame
rownames(MOA.df) <- MOA.df$Name


load('client-side/output/analyze.L1000.data.R.output/A549.analyze.L1000.RData')
cp.fit.df$fdr <- p.adjust(cp.fit.df$p.value.positive.effect,method='bonferroni')
sig.cp.fit.df <- cp.fit.df[cp.fit.df$fdr < cut.off,]
sig.cp.fit.df <- sig.cp.fit.df[order(sig.cp.fit.df$effect.size,decreasing = TRUE),]
flag <- grepl(x=sig.cp.fit.df$pert_iname,pattern = 'BRD') | sig.cp.fit.df$pert_iname %in% c(HDACi,CDKi,RAFi,HSPi,DNMTi)
sig.cp.fit.df <- sig.cp.fit.df[!flag,]
