load('client-side//output//analyze.L1000.data.R.output//A375.analyze.L1000.RData')
require(data.table)
############### comparision with clue.io query, touchstone compounds ################
data                     <- fread(input='client-side/external.data/from.LINCS//A375.lincs.query.cp.txt',header=TRUE) %>% as.data.frame
tmp                      <- ddply(data,.(Name),function(x) data.frame(score=median(x$Score)))
rownames(tmp)            <- tmp$Name
tmp                      <- tmp[order(tmp$score,decreasing = TRUE),]
A375.cp.lincs.fit.df     <- tmp


common.cp            <- intersect(A375.cp.lincs.fit.df$Name,cp.fit.df.24$cp.name) %>% as.character
A375.cp.lincs.fit.df <- A375.cp.lincs.fit.df[A375.cp.lincs.fit.df$Name %in% common.cp,]

cp.sig.df.24         <- cp.fit.df.24[cp.fit.df.24$fdr < 0.01,]
cp.sig.df.24         <- cp.sig.df.24[cp.sig.df.24$cp.name %in% common.cp,]

idx <- nrow(cp.sig.df.24)

lincs.cp <- A375.cp.lincs.fit.df$Name[1:idx]
my.cp    <- cp.sig.df.24$cp.name %>% as.character


HDACi <- read.table('client-side/external.data/from.LINCS//LINCS.HDACi.txt',header=FALSE)$V1 %>% as.character %>% unique
CDKi  <- read.table('client-side/external.data/from.LINCS/LINCS.CDKi.txt',header=FALSE)$V1 %>% as.character  %>% unique
RAFi  <- read.table('client-side/external.data/from.LINCS/LINCS.RAFi.txt',header=FALSE)$V1 %>% as.character  %>% unique
HSPi  <- read.table('client-side/external.data/from.LINCS/LINCS.HSPi.txt',header=FALSE)$V1 %>% as.character  %>% unique
DNMTi  <- read.table('client-side/external.data/from.LINCS/LINCS.DNMTi.txt',header=FALSE)$V1 %>% as.character  %>% unique


HDACi <- intersect(HDACi,common.cp)
CDKi  <- intersect(CDKi,common.cp)
RAFi  <- intersect(RAFi,common.cp)
HSPi  <- intersect(HSPi,common.cp)
DNMTi <- intersect(DNMTi,common.cp)

#intersect(lincs.cp,CDKi)
#intersect(my.cp,CDKi)
CDKi.my.p.value    <- phyper( q = intersect(my.cp,CDKi) %>% length    - 1   , m=length(CDKi), n=length(common.cp)-length(CDKi), k=length(my.cp), lower.tail = F)   
CDKi.lincs.p.value <- phyper(q = intersect(lincs.cp,CDKi) %>% length - 1  , m=length(CDKi), n=length(common.cp)-length(CDKi), k=length(my.cp), lower.tail = F)   



#intersect(my.cp,HDACi)
#intersect(lincs.cp,HDACi)
HDACi.my.p.value    <- phyper( q = intersect(my.cp,HDACi) %>% length   - 1   , m=length(HDACi), n=length(common.cp)-length(HDACi), k=length(my.cp), lower.tail = F)   
HDACi.lincs.p.value <- phyper(q = intersect(lincs.cp,HDACi) %>% length - 1  , m=length(HDACi), n=length(common.cp)-length(HDACi), k=length(my.cp), lower.tail = F)   


#intersect(my.cp,RAFi)
#intersect(lincs.cp,RAFi)
RAFi.my.p.value    <- phyper( q = intersect(my.cp,RAFi) %>% length    - 1   , m=length(RAFi), n=length(common.cp)-length(RAFi), k=length(my.cp), lower.tail = F)   
RAFi.lincs.p.value <- phyper(q = intersect(lincs.cp,RAFi) %>% length - 1  , m=length(RAFi), n=length(common.cp)-length(RAFi), k=length(my.cp), lower.tail = F)   



# intersect(my.cp,HSPi)
# intersect(lincs.cp,HSPi)
HSPi.my.p.value    <- phyper( q = intersect(my.cp,HSPi) %>% length    - 1   , m=length(HSPi), n=length(common.cp)-length(HSPi), k=length(my.cp), lower.tail = F)   
HSPi.lincs.p.value <- phyper(q = intersect(lincs.cp,HSPi) %>% length - 1  , m=length(HSPi), n=length(common.cp)-length(HSPi), k=length(my.cp), lower.tail = F)   



#intersect(my.cp,DNMTi)
#intersect(lincs.cp,DNMTi)
DNMTi.my.p.value    <- phyper( q = intersect(my.cp,DNMTi) %>% length    - 1   , m=length(DNMTi), n=length(common.cp)-length(DNMTi), k=length(my.cp), lower.tail = F)   
DNMTi.lincs.p.value <- phyper(q = intersect(lincs.cp,DNMTi) %>% length - 1  , m=length(DNMTi), n=length(common.cp)-length(DNMTi), k=length(my.cp), lower.tail = F)   

