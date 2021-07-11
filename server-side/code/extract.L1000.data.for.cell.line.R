require(cmapR)
require(dplyr)
require(data.table)
require(GSVA)


target.cell.line <- 'HELA'
# target.cell.line <- 'HEPG2'
# target.cell.line <- 'A549'
# target.cell.line <- 'MCF7'
# 
# tmp                <- fread(input='L1000.land.mark.gene.tsv.txt',header=TRUE) %>% as.data.frame
# land.mark.gene.id  <- tmp[,1] %>% as.character

land.mark.gene.id      <- system('cut -f 1 L1000.land.mark.gene.tsv.txt',wait = TRUE,intern = TRUE) %>% as.character
land.mark.gene.id      <- land.mark.gene.id[-1]

inst.info.df           <- fread(input='PhaseI/GSE92742_Broad_LINCS_inst_info.txt',header = TRUE,stringsAsFactors = FALSE) %>% as.data.frame
flag                   <- inst.info.df[,'cell_id'] == target.cell.line
if(sum(flag) == 0) {
    cell.line.inst.info.df.back <- NULL
    cell.line.inst.matrix.back  <- NULL
}else{
    inst.id                <- inst.info.df[flag,'inst_id']
    cell.line.inst.info.df <- inst.info.df[flag,]
    rs                     <- parse.gctx('PhaseI/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx',cid=inst.id,rid=land.mark.gene.id)
    cell.line.inst.matrix  <- rs@mat
    cell.line.inst.info.df.back <- cell.line.inst.info.df
    cell.line.inst.matrix.back  <- cell.line.inst.matrix
}



col.name <- c('inst_id',  'rna_plate',	'rna_well',	'pert_id',	'pert_iname',	'pert_type',	'pert_dose',	'pert_dose_unit',	'pert_time',	'pert_time_unit',	'cell_id')
inst.info.df           <- fread(input='PhaseII/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt',header = TRUE,stringsAsFactors = FALSE) %>% as.data.frame
flag                   <- inst.info.df[,'cell_id'] == target.cell.line
if(sum(flag)!=0) {
    inst.id                <- inst.info.df[flag,'inst_id']
    cell.line.inst.info.df <- inst.info.df[flag,]
    rs                     <- parse.gctx('PhaseII/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx',cid=inst.id,rid=land.mark.gene.id)
    cell.line.inst.matrix  <- rs@mat

    cell.line.inst.info.df$pert_mfc_id    <- NULL
    colnames(cell.line.inst.info.df)[3:4] <- c('rna_plate','rna_well')
    cell.line.inst.info.df <- cell.line.inst.info.df[,col.name]

    cell.line.inst.info.df <- rbind(cell.line.inst.info.df,cell.line.inst.info.df.back)
    cell.line.inst.matrix  <- cbind(cell.line.inst.matrix, cell.line.inst.matrix.back)
} else{
    cell.line.inst.info.df <- cell.line.inst.info.df.back
    cell.line.inst.matrix  <- cell.line.inst.matrix.back
}
#ssgsea.score.matrix    <- gsva(expr = cell.line.inst.matrix[land.mark.gene.id,],gset.idx.list = list(IFN.gamma = genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']],IFN.alpha = genesets[['HALLMARK_INTERFERON_ALPHA_RESPONSE']]),method = 'ssgsea',ssgsea.norm=FALSE)

save(file=paste(target.cell.line,".lincs.RData",sep=''),list = c('cell.line.inst.info.df','cell.line.inst.matrix'))






#save(file='A375.sig.RData',list = c('cell.line.sig.info.df','cell.line.sig.matrix'))
# sig.info.df            <- fread(input='PhaseI/GSE92742_Broad_LINCS_sig_info.txt',header = TRUE,stringsAsFactors = FALSE) %>% as.data.frame
# flag                   <- sig.info.df[,'cell_id'] == target.cell.line
# sig.id                 <- sig.info.df[flag,'sig_id']
# cell.line.sig.info.df  <- sig.info.df[flag,]
# rs                     <- parse.gctx('PhaseI/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx',cid=sig.id,rid=land.mark.gene.id)
# cell.line.sig.matrix   <- rs@mat
