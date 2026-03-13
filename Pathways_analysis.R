#Packages 
library(org.Mm.eg.db)
library(org.Dr.eg.db)
library(org.At.tair.db)
library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(geneslator)
library(clusterProfiler)
library(enrichplot)
library(xlsx)
library(openxlsx)
library(readxl)
#-----------------------------------------------------------------------------#
setwd("./")
#-----------------------------------------------------------------------------#

kegg_analysis <- function(df, input_col, from_type, pckg_org, pckg_gns, kegg_org){
  # Gene mapping
  if (pckg_gns == "org.Scerevisiae.db"){
    clf <- clusterProfiler::bitr(df[,input_col], fromType = from_type, toType = "ENSEMBL", OrgDb = pckg_org)
    gns <- geneslator::select(
      get(pckg_gns),
      keys = df[[input_col]],
      keytype = from_type,
      columns = "ENSEMBL")
    kk1_gns <- enrichKEGG(gene = gns$ENSEMBL, organism = kegg_org, pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1)
    kk1_clf <- enrichKEGG(gene = clf$ENSEMBL, organism = kegg_org, pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1)
  } else if(pckg_gns == "org.Athaliana.db"){
    clf <- clusterProfiler::bitr(df[,input_col], fromType = from_type, toType = "TAIR", OrgDb = pckg_org)
    gns <- geneslator::select(
      get(pckg_gns),
      keys = df[[input_col]],
      keytype = from_type,
      columns = "ENSEMBL")
    kk1_gns <- enrichKEGG(gene = gns$ENSEMBL, organism = kegg_org, pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1)
    kk1_clf <- enrichKEGG(gene = clf$TAIR, organism = kegg_org, pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1)
  } else {
    clf <- clusterProfiler::bitr(df[,input_col], fromType = from_type, toType = "ENTREZID", OrgDb = pckg_org)
    gns <- geneslator::select(
      get(pckg_gns),
      keys = df[[input_col]],
      keytype = from_type,
      columns = "ENTREZID")
    kk1_gns <- enrichKEGG(gene = gns$ENTREZID, organism = kegg_org, pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1)
    kk1_clf <- enrichKEGG(gene = clf$ENTREZID, organism = kegg_org, pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1)
  }
 
  # Enrichment KEGG

  # Dataframe KEGG
  df_kk_gns <- data.frame(
    KEGG_Description = kk1_gns@result$Description, 
    GeneRatio_gns = kk1_gns@result$GeneRatio, 
    KEGG_Pvalue_gns = kk1_gns@result$pvalue,
    GeneID_gns = kk1_gns@result$geneID
  )
  df_kk_clf <- data.frame(
    KEGG_Description = kk1_clf@result$Description, 
    GeneRatio_clf = kk1_clf@result$GeneRatio, 
    KEGG_Pvalue_clf = kk1_clf@result$pvalue,
    GeneID_clf = kk1_clf@result$geneID
  )
  # Full join
  tot_kk <- merge(df_kk_gns, df_kk_clf, by = "KEGG_Description", all = TRUE)
  # Pathways differences
  diff_pathways <- setdiff(
    tot_kk$KEGG_Description[!is.na(tot_kk$KEGG_Pvalue_gns)],
    tot_kk$KEGG_Description[!is.na(tot_kk$KEGG_Pvalue_clf)]
  )
  if(length(diff_pathways) == 0){
    # GeneRatio calculation
    tot_kk$gene_gns <- as.numeric(gsub("/.*","", tot_kk$GeneRatio_gns))
    tot_kk$gene_clf <- as.numeric(gsub("/.*","", tot_kk$GeneRatio_clf))
    tot_kk$gene_ratio_diversity <- ifelse(tot_kk$gene_gns != tot_kk$gene_clf, "YES", "NO")
    tot_kk_yes <- tot_kk[tot_kk$gene_ratio_diversity == "YES",]
    tot_kk_yes$org_mig <- ifelse(tot_kk_yes$gene_gns < tot_kk_yes$gene_clf, "YES", "NO")
    tot_kk_org <- tot_kk_yes[tot_kk_yes$org_mig == "YES",]
    tot_kk_gns <- tot_kk_yes[tot_kk_yes$org_mig == "NO",]
    tot_kk_gns$diff <- tot_kk_gns$gene_gns - tot_kk_gns$gene_clf
    diff_mean <- mean(tot_kk_gns$diff)
    # Differenze tra geni
    trova_differenze <- function(lista1, lista2){
      geni1 <- unlist(strsplit(lista1,"/"))
      geni2 <- unlist(strsplit(lista2,"/"))
      list(
        only_in_gns = paste(setdiff(geni1, geni2), collapse="/"),
        only_in_clf = paste(setdiff(geni2, geni1), collapse="/")
      )
    }
    differenze <- mapply(trova_differenze, tot_kk_org$GeneID_gns, tot_kk_org$GeneID_clf, SIMPLIFY = FALSE)
    tot_kk_org$gns_only_genes <- sapply(differenze, function(x) x$only_in_gns)
    tot_kk_org$clf_only_genes <- sapply(differenze, function(x) x$only_in_clf)
    differences_df <- tot_kk_org[tot_kk_org$gns_only_genes != "" | tot_kk_org$clf_only_genes != "", ]
    
    file_out <- paste0(
      "kegg_", kegg_org, "_org_vs_geneslator_", from_type, ".xlsx"
    )
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = from_type)
    writeData(wb, sheet = from_type, x = tot_kk_yes)
    addWorksheet(wb, sheetName = paste0(from_type, "_diff"))
    writeData(wb, sheet = paste0(from_type, "_diff"), x = differences_df)
    saveWorkbook(wb, file_out, overwrite = TRUE)
    Top_mapped_pathways <- round((nrow(tot_kk_gns)/nrow(tot_kk))*100,2)
    tot_kk$diff <- tot_kk$gene_gns - tot_kk$gene_clf
    wx_test <- wilcox.test(tot_kk$diff, mu = 0, alternative = "greater")
  } else {
    # extra Pathways
    pathways_extra <- tot_kk[tot_kk$KEGG_Description %in% diff_pathways, ]
    pathways_extra$origine <- ifelse(!is.na(pathways_extra$KEGG_Pvalue_gns), "gns", "clf")
    cat("=== Extra pathways ===\n")
    print(pathways_extra$KEGG_Description)
    write.xlsx(pathways_extra, paste0("pathways_extra_", kegg_org, ".xlsx"),
               sheetName = "extra_pathways", row.names = FALSE)
    Top_mapped_pathways <- NA
    wx_test <- list(p.value = NA)
  }
  summary <- data.frame(
    Top_mapped_pathways = Top_mapped_pathways,
    Wilcoxon_pvalue = wx_test$p.value,
    average_n_of_extra_genes <- diff_mean
  )
  return(summary)
}

#Homo sapiens
GeneslatorDb("Homo sapiens")
#SYMBOL
human<- readRDS("Human/example.symbolid.RDS")
human <- unique(human)

human_symbol <- analisi_pathways(df = human,
                                 input_col = "hgnc_symbol",
                                 from_type = "SYMBOL",
                                 pckg_org = "org.Hs.eg.db",
                                 pckg_gns = "org.Hsapiens.db",
                                 kegg_org = "hsa")


#ENSEMBL
human<- readRDS("Human/list.ensembl.RDS")
human <- unique(human)

human_ensembl <- analisi_pathways(df = human,
                                  input_col = "Ensembl.gene.ID",
                                  from_type = "ENSEMBL",
                                  pckg_org = "org.Hs.eg.db",
                                  pckg_gns = "org.Hsapiens.db",
                                  kegg_org = "hsa")
#Mus musculus
GeneslatorDb("Mus musculus")
#SYMBOL
mouse1 <- readRDS("Mouse/example.symbolid.RDS") 
mouse1 <- unique(mouse1)
rownames(mouse1) <- NULL
mouse1 <- as.data.frame(mouse1[-29759,])
rownames(mouse1) <- NULL
colnames(mouse1) <- "gene_name"

mouse_symbol <- kegg_analysis(df = mouse1,
                                input_col = "gene_name",
                                from_type = "SYMBOL",
                                pckg_org = "org.Mm.eg.db",
                                pckg_gns = "org.Mmusculus.db",
                                kegg_org = "mmu")

#ENSEMBL
mouse1 <- readRDS("Mouse/example.ensembl.RDS")
mouse1$Ensembl <- as.character(mouse1$Ensembl)

mouse_ensembl <- kegg_analysis(df = mouse1,
                                    input_col = "Ensembl",
                                    from_type = "ENSEMBL",
                                    pckg_org = "org.Mm.eg.db",
                                    pckg_gns = "org.Mmusculus.db",
                                    kegg_org = "mmu")

#Danio rerio
GeneslatorDb("Danio rerio")
#SYMBOL
zebr <- read.csv("Zebrafish/conte_zebrafish_symbol.txt", sep=",", header=FALSE) 
colnames(zebr)[1] <- "SYMBOL"
zebr <- as.data.frame(zebr$SYMBOL)
colnames(zebr)[1] <- "SYMBOL"

zebrafish_symbol <- kegg_analysis(df = zebr,
                                input_col = "SYMBOL",
                                from_type = "SYMBOL",
                                pckg_org = "org.Dr.eg.db",
                                pckg_gns = "org.Drerio.db",
                                kegg_org = "dre")

#ENSEMBL
zebr <- read_excel("Zebrafish/GSE274820_Repository_AnalysisData_2019_A5806_DESeq.xlsx", sheet = 3) 
colnames(zebr)[1] <- "ENSEMBL"
zebr <- as.data.frame(zebr$ENSEMBL)
colnames(zebr)[1] <- "ENSEMBL"

zebrafish_ensembl <- kegg_analysis(df = zebr,
                                 input_col = "ENSEMBL",
                                 from_type = "ENSEMBL",
                                 pckg_org = "org.Dr.eg.db",
                                 pckg_gns = "org.Drerio.db",
                                 kegg_org = "dre")


#Saccharomyces cerevisiae
GeneslatorDb("Saccharomyces cerevisiae")
#ENTREZID
scer1 <- read.csv("Yeast/Saccharomyces_cerevisiae_ENTREZID_GPL90.txt", sep = "\t", header = FALSE) 
colnames(scer1) <- "GeneID"
scer1$GeneID <- as.character(scer1$GeneID)

yeast_entrezid <- kegg_analysis(df = scer1,
                                 input_col = "GeneID",
                                 from_type = "ENTREZID",
                                 pckg_org = "org.Sc.sgd.db",
                                 pckg_gns = "org.Scerevisiae.db",
                                 kegg_org = "sce")

#Rattus norvegicus
GeneslatorDb("Rattus norvegicus")
#SYMBOL
rn <- readRDS("Rattus/example.symbolid.RDS")
rn <- unique(rn)

rattus_symbol <- kegg_analysis(df = rn,
                              input_col = "gene",
                              from_type = "SYMBOL",
                              pckg_org = "org.Rn.eg.db",
                              pckg_gns = "org.Rnorvegicus.db",
                              kegg_org = "rno")

#ENSEMBL
rn <- readRDS("Rattus/ensembl.RDS")
rn <- unique(rn)

rattus_ensembl <- kegg_analysis(df = rn,
                               input_col = "Ensembl.gene.ID",
                               from_type = "ENSEMBL",
                               pckg_org = "org.Rn.eg.db",
                               pckg_gns = "org.Rnorvegicus.db",
                               kegg_org = "rno")

#Arabidopsis thaliana
GeneslatorDb("Arabidopsis thaliana")
#SYMBOL
ara1 <- read.csv("Arabidopsis/GSM8858406_TPM_Ler_R1.csv", sep=",") #pulire gli AC, NC, RA
ara1 <- as.data.frame(ara1$Ara_gene)
ara1$`ara1$Ara_gene` <- gsub("\\|.*", "", ara1$`ara1$Ara_gene`)
colnames(ara1) <- "Gene"

arabidopsis_symbol <- kegg_analysis(df = ara1,
                                       input_col = "Gene",
                                       from_type = "SYMBOL",
                                       pckg_org = "org.At.tair.db",
                                       pckg_gns = "org.Athaliana.db",
                                       kegg_org = "ath")

#ENTREZID
ara1 <- read.csv("Arabidopsis/Arabidopsis_thaliana_EntrezID_GPL198.txt", sep="\t", header= FALSE) 
colnames(ara1) <- "GeneID"
ara1$GeneID <- as.character(ara1$GeneID)

arabidopsis_entrezid <- kegg_analysis(df = ara1,
                                    input_col = "GeneID",
                                    from_type = "ENTREZID",
                                    pckg_org = "org.At.tair.db",
                                    pckg_gns = "org.Athaliana.db",
                                    kegg_org = "ath")

