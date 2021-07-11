require(GSA)
require(plyr)
require(dplyr)
require(limma)
require(hgu133a.db)

####### L1000 gene entrez id
load('server-side//RData//lincs/lincs_signatures_cmpd.RData')
L1000.gene.entrez.id <- rownames(lincs_signatures_cmpd) %>% as.character


####### Stable gene ids and their mapping between ensemble and entrez
hgnc_complete_set  <- read.delim("~/Project/CancerImmunoTherapy-TCGA/client-side/external.data/hgnc.annotation/hgnc_complete_set.txt", stringsAsFactors=FALSE)

df1                 <- table(hgnc_complete_set$entrez_id)       %>% as.data.frame %>% arrange(desc(Freq)) #  entrez ids mapping to multiple ensemble id
exclude.entrez.id   <- df1$Var1[df1$Freq > 1] %>% as.character
exclude.entrez.id   <- c(exclude.entrez.id,'')
df2                 <- table(hgnc_complete_set$ensembl_gene_id) %>% as.data.frame %>% arrange(desc(Freq)) #  ensemble ids mapping to multiple entrez id
exclude.ensemble.id <- df2$Var1[df2$Freq > 1] %>% as.character

hgnc_complete_set   <- hgnc_complete_set[,c('ensembl_gene_id','entrez_id','locus_group')]
flag                <- (hgnc_complete_set$ensembl_gene_id %in% exclude.ensemble.id) | (hgnc_complete_set$entrez_id %in% exclude.entrez.id)
hgnc_complete_set   <- hgnc_complete_set[!flag,]
hgnc_complete_set   <- hgnc_complete_set[complete.cases(hgnc_complete_set),]

ensemble2entrez        <- hgnc_complete_set$entrez_id
names(ensemble2entrez) <- hgnc_complete_set$ensembl_gene_id

protein.coding.gene.ensemble.id   <- hgnc_complete_set$ensembl_gene_id[hgnc_complete_set$locus_group == 'protein-coding gene'] %>% unique
protein.coding.gene.entrez.id     <- hgnc_complete_set$entrez_id[hgnc_complete_set$locus_group == 'protein-coding gene'] %>% unique




####### msigdb genesets

msigdb          <-  GSA.read.gmt("client-side/external.data/MSigDB/h.all.v6.1.entrez.gmt")
genesets        <-  msigdb$genesets
names(genesets) <-  msigdb$geneset.names 



################ maaping GPL96 probe_set id to entrez gene id



process.entrez.gene.id <- function(x) {
  if(length(x) == 1){
    if(x == '' | is.na(x)){
      return(NA)
    } else {
      return(x)  
    }    
  }
  
  if(sum(x %in% protein.coding.gene.entrez.id) == 1){ # well, assign to protein coding genes since they are higly expressed
    return(x[x %in% protein.coding.gene.entrez.id]) 
  }else{
    return(NA)
  }
}

split.string <- function(x){
  x   <- gsub(pattern = ' ',replacement = '',x = x)
  tmp <- strsplit(x=x,split = '///')  %>% unlist 
  if(length(tmp) == 0) {
    tmp <- NA
  }
  tmp
}


##### for GPL96 platform, hgu133A

l1                   <- sapply(as.list(hgu133aENTREZID),process.entrez.gene.id) # hgu133A package annotation 
GPL96.probe.2.entrez <- l1[is.na(l1) == FALSE]

tmp                  <- read.delim("~/Project/panCancer.Metastatic.microenvironment/client-side/meta.data/GPL96.mart_export.txt", dec=",", stringsAsFactors=FALSE) # biomart annotation
tmp                  <- tmp[complete.cases(tmp),]
tmp$ENTREZ_GENE_ID   <- as.character(tmp$ENTREZ_GENE_ID)
df                   <- ddply(tmp,.(ID),function(x) process.entrez.gene.id(x$ENTREZ_GENE_ID))
df                   <- df[complete.cases(df),]
df                   <- df[(df$ID %in% names(GPL96.probe.2.entrez) == FALSE),] # impute missed probesets
l1                   <- df$V1
names(l1)            <- df$ID
GPL96.probe.2.entrez <- c(GPL96.probe.2.entrez,l1)

tmp                <- read.delim("~/Project/panCancer.Metastatic.microenvironment/client-side/meta.data/GPL96.probe.2.entrez.txt", dec=",", stringsAsFactors=FALSE) # GEO annotation
tmp                <- tmp[,c('ID','ENTREZ_GENE_ID')]
rs                 <- lapply(tmp$ENTREZ_GENE_ID,split.string)
l1                 <- sapply(rs,process.entrez.gene.id)
names(l1)          <- tmp$ID
l1                 <- l1[is.na(l1) == FALSE]
l1                 <- l1[(names(l1) %in% names(GPL96.probe.2.entrez)) == FALSE] # impute missed probesets
GPL96.probe.2.entrez <- c(GPL96.probe.2.entrez,l1)

#### save data
save(file='client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData',list=c('L1000.gene.entrez.id','hgnc_complete_set','ensemble2entrez','genesets','protein.coding.gene.ensemble.id','protein.coding.gene.entrez.id','GPL96.probe.2.entrez'))

# gene2ensembl   <- read.delim("~/Project/CancerImmunoTherapy-TCGA/client-side/meta.data/gene2ensembl", stringsAsFactors=FALSE)
# flag           <- grepl(x=gene2ensembl$RNA_nucleotide_accession.version,pattern = 'NM') | grepl(x=gene2ensembl$RNA_nucleotide_accession.version,pattern = 'XM')
# gene2ensembl   <- gene2ensembl[flag,]
# data           <- unique(gene2ensembl[gene2ensembl$tax_id == 9606,c('Ensembl_gene_identifier','GeneID')])
# 
# df1 <- ddply(data, .(Ensembl_gene_identifier),function(x) unique(x$GeneID)                  %>% length)
# df2 <- ddply(data, .(GeneID),                 function(x) unique(x$Ensembl_gene_identifier) %>% length)
# 
# exclude.gene.id <- df1$Ensembl_gene_identifier[df1$V1 > 1]
# exclude.ref.id  <- df2$GeneID[df2$V1 > 1]
# 
# 
# data           <- data[((data$Ensembl_gene_identifier %in% exclude.gene.id) | (data$GeneID %in% exclude.ref.id)) == FALSE ,]
# ensemble2entrez        <- data$GeneID
# ensemble2entrez        <- as.character(ensemble2entrez)
# names(ensemble2entrez) <- data$Ensembl_gene_identifier
# ensemble2entrez        <- c(ensemble2entrez,ENSG00000127603 = '23499')
# require(org.Hs.eg.db)
# ensemble.to.entrez.mapping       <- revmap(org.Hs.egENSEMBL) %>% as.list
# cnt                              <- sapply(ensemble.to.entrez.mapping,length)
# ensemble.to.entrez.mapping       <- ensemble.to.entrez.mapping[cnt == 1]
# x                                <- unlist(ensemble.to.entrez.mapping)
# x                                <- x[ (names(x) %in% names(ensemble2entrez)) == FALSE]
# ensemble2entrez                  <- c(ensemble2entrez,x)


# data                     <- read.csv("~/Project/CancerImmunoTherapy-TCGA/client-side/meta.data/TCGA.tumor.purity.csv", stringsAsFactors=FALSE)
# data$Sample.ID           <- gsub(x=data$Sample.ID ,pattern='A$',replacement = '')
# rownames(data)           <- data$Sample.ID
# TCGA.tumor.purity        <- data$ESTIMATE
# names(TCGA.tumor.purity) <- rownames(data)
# 



# # here I ONLY consider protein coding genes, their IDs are much more stable!
# 
# hgnc_complete_set                 <- read.delim("~/Project/CancerImmunoTherapy-TCGA/client-side/meta.data/hgnc_complete_set.txt", stringsAsFactors=FALSE)
# protein.coding.gene.ensemble.id   <- hgnc_complete_set$ensembl_gene_id[hgnc_complete_set$locus_group == 'protein-coding gene'] %>% unique
# protein.coding.gene.ensemble.id   <- protein.coding.gene.ensemble.id[protein.coding.gene.ensemble.id != '']
# 
# hgnc_complete_set                 <- hgnc_complete_set[,c('ensembl_gene_id','entrez_id')]
# hgnc_complete_set                 <- hgnc_complete_set[hgnc_complete_set$ensembl_gene_id != '',]
# hgnc_complete_set$ensembl_gene_id <- as.character(hgnc_complete_set$ensembl_gene_id)
# hgnc_complete_set$entrez_id       <- as.character(hgnc_complete_set$entrez_id)
# 
# df                     <- table(hgnc_complete_set$ensembl_gene_id) %>% as.data.frame %>% arrange(desc(Freq))
# exclude.gene           <- df$Var1[df$Freq > 1] # exclude genes mapping to multiple entrez gene ids
# hgnc_complete_set      <- hgnc_complete_set[hgnc_complete_set$ensembl_gene_id!=exclude.gene,]
# ensemble2entrez        <- hgnc_complete_set$entrez_id
# names(ensemble2entrez) <- hgnc_complete_set$ensembl_gene_id

