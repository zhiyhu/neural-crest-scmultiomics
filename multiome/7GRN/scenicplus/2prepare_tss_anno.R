# Prepare TSS annotation 
# Zhiyuan Hu
# 21 Jan 2023
# last modified 21 Jan 2023
wkdir="multiome/analysis_newref/"
seu_var <- read.delim(paste0(wkdir,"preprocessing/rds/seu_featuredata.tsv"))

pybiomart <- read.csv("multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105.csv",
                      row.names = 1)

seu_var$gene_id[!seu_var$gene_id%in% pybiomart$Gene.stable.ID]
# [1] "foxd3-citrine" "foxd3-mCherry"

pybiomart <- pybiomart[pybiomart$Transcript.type == 'protein_coding',]
head(pybiomart)

idx <- !grepl("CHR",pybiomart$Chromosome.scaffold.name) & !grepl("GL",pybiomart$Chromosome.scaffold.name) & !grepl("JH",pybiomart$Chromosome.scaffold.name)& 
  !grepl("MT",pybiomart$Chromosome.scaffold.name) & !grepl("KN",pybiomart$Chromosome.scaffold.name)& !grepl("KZ",pybiomart$Chromosome.scaffold.name)
pybiomart <- pybiomart[idx,]

pybiomart <- pybiomart[pybiomart$Gene.stable.ID %in% seu_var$gene_id, ]
pybiomart$seu_gene_name <- seu_var$gene_name_unique[match(pybiomart$Gene.stable.ID, seu_var$gene_id)]

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
write.csv(pybiomart, "multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105_geneNameMatched.csv",row.names = FALSE)


