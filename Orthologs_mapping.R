#Ortologhs

#Packages
library(geneslator)
library(orthogene)
library(homologene)
library(biomaRt)
library(babelgene)
library(gprofiler2)
library(readxl)

#------------------------------------------------------------#
setwd("./")
#------------------------------------------------------------#
#Function
orthologs_mapping <- 
  function(df, pckg, col_from, col_input, specie, taxa, 
           species_name, col_from_lowercase, species_code){
    #geneslator
    annot <- geneslator::select(pckg,
                                keys = df[,col_input],
                                columns = c(col_from, "ORTHOHUMAN"),
                                keytype = col_from)
    total_geneslator <- annot[!is.na(annot$ORTHOHUMAN),]
    total_geneslator <- total_geneslator[!duplicated(total_geneslator[,col_from]),]
    geneslator <- round(nrow(total_geneslator)/nrow(df)*100,2)
    #orthogene - gprofiler
    total_orthogene_gprofiler <- orthogene::convert_orthologs(gene_df = df,
                                                              gene_input = col_input, 
                                                              gene_output = "columns", 
                                                              input_species = specie,
                                                              output_species = "human",
                                                              non121_strategy = "keep_both_species", 
                                                              method = "gprofiler")
    total_orthogene_gprofiler <- total_orthogene_gprofiler[!duplicated(total_orthogene_gprofiler$input_gene),]
    total_orthogene_gprofiler <- total_orthogene_gprofiler[!is.na(total_orthogene_gprofiler$ortholog_gene),]
    ort_gp <- round(nrow(total_orthogene_gprofiler)/nrow(df)*100,2)
    #orthogene - babelgene
    total_orthogene_babelgene <- orthogene::convert_orthologs(gene_df = df,
                                                              gene_input = col_input, 
                                                              gene_output = "columns", 
                                                              input_species = specie,
                                                              output_species = "human",
                                                              non121_strategy = "keep_both_species", 
                                                              method = "babelgene")
    total_orthogene_babelgene <- total_orthogene_babelgene[!duplicated(total_orthogene_babelgene$input_gene),]
    total_orthogene_babelgene <- total_orthogene_babelgene[!is.na(total_orthogene_babelgene$ortholog_gene),]
    ort_bab <- round(nrow(total_orthogene_babelgene)/nrow(df)*100,2)
    #Homologene
    human_orthologs <- homologene(df[,col_input], inTax = taxa, outTax = 9606)
    human_orthologs <- human_orthologs[,c(1,2)]
    human_orthologs <- human_orthologs[!is.na(human_orthologs$`9606`),]
    human_orthologs <- human_orthologs[!duplicated(human_orthologs[,1]),]
    homologene <- round(nrow(human_orthologs)/nrow(df)*100,2)
    #babelgene
    if(species_name == "Rattus norvegicus" & col_from == "ENSEMBL"){
      babelgene <- 0
    } else{
      orthologs <- orthologs(genes = df[,col_input], 
                           species = species_name, human=FALSE)
    orthologs <- orthologs[,c(col_from_lowercase, paste0("human_", col_from_lowercase))]
    orthologs <- orthologs[!is.na(orthologs[,paste0("human_", col_from_lowercase)]),]
    orthologs <- orthologs[!duplicated(orthologs[,col_from_lowercase]),]
    babelgene <- round(nrow(orthologs)/nrow(df)*100,2)
    }
   
    #gprofiler2
    if(col_from == "ENTREZID"){
      if(specie == "drosophila"){
        result <- gorth(query = df[,col_input],
                        source_organism = species_code,
                        target_organism = "hsapiens", numeric_ns = "ENTREZGENE")
      } else {
      result <- gorth(query = df[,col_input],
                      source_organism = species_code,
                      target_organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC")
    }} else {
      result <- gorth(query = df[,col_input],
                      source_organism = species_code,
                      target_organism = "hsapiens")
    }
    result <- result[,c("input", "ortholog_name")] 
    result <- result[!is.na(result$ortholog_name),]
    result <- result[!duplicated(result$input),]
    gprofiler <- round(nrow(result)/nrow(df)*100,2)
    #total
    result<- data.frame("Specie"=species_name,
                        "geneslator"=geneslator, 
                        "orthogene_gprofiler"=ort_gp,
                        "orthogene_babelgene"=ort_bab, 
                        "homologene"=homologene, 
                        "gprofiler2"=gprofiler,
                         "babelgene" = babelgene)
    return(result)
  }

#----------------------------------------------------------------------------------------#
#SYMBOL
#Mouse
GeneslatorDb("Mus musculus")
mouse1 <- readRDS("mouse/example.symbolid.RDS") 
mouse1 <- unique(mouse1)
rownames(mouse1) <- NULL
mouse1 <- as.data.frame(mouse1[-29759,])
rownames(mouse1) <- NULL
colnames(mouse1) <- "gene_name"
result_mouse <- orthologs_mapping(df = mouse1,
                                      pckg = org.Mmusculus.db,
                                      col_from= "SYMBOL",
                                      col_input= "gene_name",
                                      specie = "mouse",
                                      taxa = 10090,
                                      species_name = "Mus musculus",
                                      col_from_lowercase = "symbol",
                                      species_code = "mmusculus")

#Zebrafish
GeneslatorDb("Danio rerio")
zebr <- read.csv("zebrafish/conte_zebrafish_symbol.txt", sep=",", header=FALSE) 
colnames(zebr)[1] <- "SYMBOL"
zebr <- as.data.frame(zebr$SYMBOL)
colnames(zebr)[1] <- "SYMBOL"
result_zebrafish <- orthologs_mapping(df = zebr,
                                      pckg = org.Drerio.db,
                                      col_from= "SYMBOL",
                                      col_input= "SYMBOL",
                                      specie = "zebrafish",
                                      taxa = 7955,
                                      species_name = "Danio rerio",
                                      col_from_lowercase = "symbol",
                                      species_code = "drerio")

#Yeast
GeneslatorDb("Saccharomyces cerevisiae")
yeast1 <- read.csv("yeast/lista_geni_yeast.txt", 
                   sep = ";", header=FALSE)
yeast1$V1 <- gsub("_", "", yeast1$V1)
yeast1 <- as.data.frame(yeast1[,c(1)])
yeast1 <- unique(yeast1)
colnames(yeast1) <- "Gene"
result_yeast <- orthologs_mapping(df = yeast1,
                                          pckg = org.Scerevisiae.db,
                                          col_from= "SYMBOL",
                                          col_input= "Gene",
                                          specie = "yeast",
                                          taxa = 4932,
                                          species_name = "Saccharomyces cerevisiae",
                                          col_from_lowercase = "symbol",
                                          species_code = "scerevisiae")

#Rat
GeneslatorDb("Rattus norvegicus")
rattus <- readRDS("Rattus/example.symbolid.RDS") 
rattus <- unique(rattus)
colnames(rattus) <- "SYMBOL"

result_rat <- orthologs_mapping(df = rattus,
                                  pckg = org.Rnorvegicus.db,
                                  col_from= "SYMBOL",
                                  col_input= "SYMBOL",
                                  specie = "rat",
                                  taxa = 10116,
                                  species_name = "Rattus norvegicus",
                                  col_from_lowercase = "symbol",
                                  species_code = "rnorvegicus")

#Drosophila
GeneslatorDb("Drosophila melanogaster")
dros <- readRDS("Drosophila/example.symbolid.RDS")
dros <- unique(dros)

result_drosophila <- orthologs_mapping(df = dros,
                                       pckg = org.Dmelanogaster.db,
                                       col_from="SYMBOL",
                                       col_input= "gene",
                                       specie = "drosophila",
                                       taxa = 7227,
                                       species_name = "Drosophila melanogaster",
                                       col_from_lowercase = "symbol",
                                       species_code = "dmelanogaster")
#Worm
GeneslatorDb("Caenorhabditis elegans")
input_df <- readRDS("Celegans/example.symbolid.RDS")
colnames(input_df)[1] <- "SYMBOL"
input_df <- unique(input_df)
result_celegans <- orthologs_mapping(df = input_df,
                                     pckg = org.Celegans.db,
                                     col_from="SYMBOL",
                                     col_input= "SYMBOL",
                                     specie = "Caenorhabditis",
                                     taxa = 6239,
                                     species_name = "Caenorhabditis elegans",
                                     col_from_lowercase = "symbol",
                                     species_code =  "celegans")
#----------------------------------------------------------------------------------------------#
#ENSEMBL
#Mouse
mouse1 <- readRDS("mouse/example.ensembl.RDS")
mouse1$Ensembl <- as.character(mouse1$Ensembl)
result_mouse <- orthologs_mapping(df = mouse1,
                                      pckg = org.Mmusculus.db,
                                      col_from= "ENSEMBL",
                                      col_input= "Ensembl",
                                      specie = "mouse",
                                      taxa = 10090,
                                      species_name = "Mus musculus",
                                      col_from_lowercase = "ensembl",
                                      species_code = "mmusculus")
#Zebrafish
zebr <- read_excel("zebrafish/GSE274820_Repository_AnalysisData_2019_A5806_DESeq.xlsx", sheet = 3) 
colnames(zebr)[1] <- "ENSEMBL"
zebr <- as.data.frame(zebr$ENSEMBL)
colnames(zebr)[1] <- "ENSEMBL"
result_zebrafish <- orthologs_mapping(df = zebr,
                                          pckg = org.Drerio.db,
                                          col_from= "ENSEMBL",
                                          col_input= "ENSEMBL",
                                          specie = "zebrafish",
                                          taxa = 7955,
                                          species_name = "Danio rerio",
                                          col_from_lowercase = "ensembl",
                                          species_code = "drerio")


#Yeast
yeast1 <- read.csv("yeast/GSE307461_raw_gene_counts.tsv", sep="\t") #pulire gli AC, NC, RA
yeast1 <- as.data.frame(yeast1$Geneid)
colnames(yeast1) <- "Geneid"
result_yeast <- orthologs_mapping(df = yeast1,
                                      pckg = org.Scerevisiae.db,
                                      col_from= "SYMBOL",
                                      col_input= "Geneid",
                                      specie = "yeast",
                                      taxa = 4932,
                                      species_name = "Saccharomyces cerevisiae",
                                      col_from_lowercase = "ensembl",
                                      species_code = "scerevisiae")
#Rat
rattus <- readRDS("Rattus/ensembl.RDS") 
rattus <- unique(rattus)
colnames(rattus) <- "ENSEMBL"

result_rat <- orthologs_mapping(df = rattus,
                                pckg = org.Rnorvegicus.db,
                                col_from= "ENSEMBL",
                                col_input= "ENSEMBL",
                                specie = "rat",
                                taxa = 10116,
                                species_name = "Rattus norvegicus",
                                col_from_lowercase = "ensembl",
                                species_code = "rnorvegicus")

#Drosophila
dros <- readRDS("Drosophila/ensembl.RDS")
dros <- unique(dros)
colnames(dros) <- "ENSEMBL"

result_drosophila <- orthologs_mapping(df = dros,
                                       pckg = org.Dmelanogaster.db,
                                       col_from="ENSEMBL",
                                       col_input= "ENSEMBL",
                                       specie = "drosophila",
                                       taxa = 7227,
                                       species_name = "Drosophila melanogaster",
                                       col_from_lowercase = "ensembl",
                                       species_code = "dmelanogaster")

#Worm
input_df <- readRDS("Celegans/ensembl.RDS")
colnames(input_df)[1] <- "ENSEMBL"
input_df <- unique(input_df)
result_celegans <- orthologs_mapping(df = input_df,
                                         pckg = org.Celegans.db,
                                         col_from="ENSEMBL",
                                         col_input= "ENSEMBL",
                                         specie = "Caenorhabditis",
                                         taxa = 6239,
                                         species_name = "Caenorhabditis elegans",
                                         col_from_lowercase = "ensembl",
                                         species_code =  "celegans")

#-----------------------------------------------------------------------------#
#ENTREZID
mouse1 <- readRDS("mouse/example.entrezid.RDS")
mouse1 <- unique(mouse1)
rownames(mouse1) <- NULL
mouse1$GeneID <- as.character(mouse1$GeneID)

result_mouse <- orthologs_mapping(df = mouse1,
                                      pckg = org.Mmusculus.db,
                                      col_from= "ENTREZID",
                                      col_input= "GeneID",
                                      specie = "mouse",
                                      taxa = 10090,
                                      species_name = "Mus musculus",
                                      col_from_lowercase = "entrez",
                                      species_code = "mmusculus")
#Zebrafish
zeb1 <- read.csv("zebrafish/Danio_rerio_EntrezID_GPL18967.txt", sep="\t") 
colnames(zeb1) <- "gene_name"
zeb1$gene_name <- as.character(zeb1$gene_name)
result_zebrafish <- orthologs_mapping(df = zeb1,
                                          pckg = org.Drerio.db,
                                          col_from= "ENTREZID",
                                          col_input= "gene_name",
                                          specie = "zebrafish",
                                          taxa = 7955,
                                          species_name = "Danio rerio",
                                          col_from_lowercase = "entrez",
                                          species_code = "drerio")
#Yeast
scer1 <- read.csv("yeast/Saccharomyces_cerevisiae_ENTREZID_GPL90.txt", sep = "\t", header = FALSE) 
colnames(scer1) <- "GeneID"
scer1$GeneID <- as.character(scer1$GeneID)
result_yeast <- orthologs_mapping(df = scer1,
                                      pckg = org.Scerevisiae.db,
                                      col_from= "ENTREZID",
                                      col_input= "GeneID",
                                      specie = "yeast",
                                      taxa = 4932,
                                      species_name = "Saccharomyces cerevisiae",
                                      col_from_lowercase = "entrez",
                                      species_code = "scerevisiae")


#Rat
rattus <- readRDS("Rattus/entrezid.RDS") 
rattus <- unique(rattus)
colnames(rattus) <- "ENTREZID"

result_rat <- orthologs_mapping(df = rattus,
                                pckg = org.Rnorvegicus.db,
                                col_from= "ENTREZID",
                                col_input= "ENTREZID",
                                specie = "rat",
                                taxa = 10116,
                                species_name = "Rattus norvegicus",
                                col_from_lowercase = "entrez",
                                species_code = "rnorvegicus")

#Drosophila
entrezdm <- readRDS("Drosophila/entrez.RDS")
entrezdm$entrez <- as.numeric(entrezdm$entrez)
entrezdm <- unique(entrezdm)

result_drosophila <- orthologs_mapping(df = entrezdm,
                                       pckg = org.Dmelanogaster.db,
                                       col_from="ENTREZID",
                                       col_input= "entrez",
                                       specie = "drosophila",
                                       taxa = 7227,
                                       species_name = "Drosophila melanogaster",
                                       col_from_lowercase = "entrez",
                                       species_code = "dmelanogaster")

#Worm
input_df <- readRDS("Celegans/entrezid.RDS")
colnames(input_df)[1] <- "ENTREZID"
input_df <- unique(input_df)

result_celegans <- orthologs_mapping(df = input_df,
                                     pckg = org.Celegans.db,
                                     col_from="ENTREZID",
                                     col_input= "ENTREZID",
                                     specie = "Caenorhabditis",
                                     taxa = 6239,
                                     species_name = "Caenorhabditis elegans",
                                     col_from_lowercase = "entrez",
                                     species_code =  "celegans")