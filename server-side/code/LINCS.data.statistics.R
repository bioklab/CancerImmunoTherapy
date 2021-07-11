require(plyr)
require(dplyr)

df1 <- fread(input= 'PhaseI/GSE92742_Broad_LINCS_inst_info.txt',header=TRUE) %>% as.data.frame
df2 <- fread(input= 'PhaseII/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt',header=TRUE) %>% as.data.frame
col.name <- c('cell_id','pert_type')
df <- rbind(df1[,col.name],df2[,col.name])
df <- table(df) %>% as.data.frame %>% print
df <- dcast(df,cell_id ~ pert_type,value.var='Freq')
write.csv(x=df,file='LINCS.data.statistics.csv',quote=FALSE,row.names=FALSE)