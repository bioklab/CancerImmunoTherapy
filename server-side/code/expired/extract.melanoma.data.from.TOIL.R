load('/home/ubuntu/chenlab_v2/chenlab_data/bulk_rna_seq/by_source/TOIL.RData')
flag          <- TOIL.sample.meta$detailed_category == 'Skin Cutaneous Melanoma'
melanoma.meta <- TOIL.sample.meta[flag,]
flag          <- melanoma.meta$sample_type != 'Solid Tissue Normal'
melanoma.meta <- melanoma.meta[flag,]

sample.meta      <- melanoma.meta
melanoma.sample  <- match(sample.meta$sample.id,colnames(TOIL.log2.fpkm.matrix))
log2.fpkm.matrix <- TOIL.log2.fpkm.matrix[,melanoma.sample]

save(file='RData/Skin Cutaneous Melanoma.primary.and.metastatic.RData',list=c('sample.meta','log2.fpkm.matrix'))