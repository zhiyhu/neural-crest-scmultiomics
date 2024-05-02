# Prepare TSS annotation 
# Zhiyuan Hu
# 21 Jan 2023
# last modified 21 Jan 2023
wkdir="/t1-data/project/tsslab/zhu/multiome/analysis_newref/"
seu_var <- read.delim(paste0(wkdir,"preprocessing/rds/seu_featuredata.tsv"))

pybiomart <- read.csv("/t1-data/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105.csv",
                      row.names = 1)

seu_var$gene_id[!seu_var$gene_id%in% pybiomart$Gene.stable.ID]
# [1] "foxd3-citrine" "foxd3-mCherry"

pybiomart <- pybiomart[pybiomart$Transcript.type == 'protein_coding',]
head(pybiomart)

idx <- !grepl("CHR",pybiomart$Chromosome.scaffold.name) & !grepl("GL",pybiomart$Chromosome.scaffold.name) & !grepl("JH",pybiomart$Chromosome.scaffold.name)& 
  !grepl("MT",pybiomart$Chromosome.scaffold.name) & !grepl("KN",pybiomart$Chromosome.scaffold.name)& !grepl("KZ",pybiomart$Chromosome.scaffold.name)
pybiomart <- pybiomart[idx,]
table(pybiomart$Chromosome.scaffold.name)

table(seu_var$gene_id %in% pybiomart$Gene.stable.ID)
# FALSE  TRUE 
# 2 27597 
table(pybiomart$Gene.stable.ID %in% seu_var$gene_id)
# FALSE  TRUE 
# 10 41653 

pybiomart <- pybiomart[pybiomart$Gene.stable.ID %in% seu_var$gene_id, ]
pybiomart$seu_gene_name <- seu_var$gene_name_unique[match(pybiomart$Gene.stable.ID, seu_var$gene_id)]
table(is.na(pybiomart$seu_gene_name))
# FALSE  TRUE 
# 50724 10913 

table(pybiomart$seu_gene_name != pybiomart$Gene.name)
# FALSE 
# 41653 

head(pybiomart)

table(pybiomart$Gene.name!=pybiomart$seu_gene_name)
# FALSE  TRUE 
# 39435  2218 
pybiomart <- pybiomart[,c("Chromosome.scaffold.name","Transcription.start.site..TSS.","Strand","seu_gene_name","Transcript.type")]

colnames(pybiomart) <- c('Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type')

which(pybiomart$Gene == "foxd3")
pybiomart[which(pybiomart$Gene == "foxd3"),]

### add the same TSS to foxd3 citrine and foxd3 mCherry
tmp = data.frame(Chromosome = c(6,6),
                Start = c(32093830,32093830),
                Strand = c(-1,-1),
                Gene = c("foxd3-citrine","foxd3-mCherry"),
                Transcript_type = c("protein_coding","protein_coding"))

pybiomart <- rbind(pybiomart, tmp)
tail(pybiomart)
write.csv(pybiomart, "/t1-data/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105_geneNameMatched.csv",row.names = FALSE)


