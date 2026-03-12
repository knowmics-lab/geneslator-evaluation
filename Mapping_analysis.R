#Packages
library(dplyr)
library(biomaRt)
library(mygene)
library(gprofiler2)
library(geneslator)
library(openxlsx)
library(org.At.tair.db)
library(org.Dr.eg.db)
library(org.Mm.eg.db)
library(org.Sc.sgd.db)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)
library(readxl)

# ============== CHECK CONVERSION FUNCTION ==============
check_conversion <- function(example, result, col_example, 
                             col_result_from, col_result_to) {
  total <- nrow(example)
  colnames(example) <- col_result_from
  not_found <- setdiff(example[[col_result_from]], result[[col_result_from]])
  not_found <- as.data.frame(not_found)
  if(dim(not_found)[1] > 0) {
    colnames(not_found) <- col_result_from
    not_found[,col_result_to] <- NA
  }
  if (nrow(not_found) == 0) not_found <- numeric(0)
  na_to <- result[is.na(result[,col_result_to]),]
  na_to_void <- result[result[,col_result_to] == "",]
  na_to <- rbind(na_to, na_to_void)
  not_na <- result[!is.na(result[,col_result_to]),]
  not_na <- not_na[not_na[,col_result_to] != "",]
  na_to <- setdiff(na_to[,col_result_from], not_na[,col_result_from])
  na_to <- as.data.frame(na_to)
  if(dim(na_to)[1] > 0) {
    colnames(na_to) <- col_result_from
    na_to[,col_result_to] <- NA
  }
  unique_df <- not_na %>%
    filter(!duplicated(not_na[,col_result_from]) & !duplicated(not_na[,col_result_from], fromLast = TRUE))
  if (nrow(na_to) == 0) na_to <- numeric(0)
  tot <- rbind(not_found, na_to)
  if(length(tot) > 0){
    na_confermati <- merge(tot, example, by=col_result_from)
    missing <- nrow(na_confermati)
  } else {
    missing <- 0
  }
  intersection_with_df <- merge(unique_df, example)
  errori <- setdiff(unique_df, intersection_with_df)
  unique_df <- setdiff(unique_df, errori)
  unique_ids <- nrow(unique_df)
  not_na1 <- merge(not_na, example)
  duplicate <- setdiff(not_na1, unique_df)
  duplicate <- duplicate[duplicate[,col_result_to] != "",]
  duplicated_genes <- length(unique(duplicate[,col_result_from]))
  summary <- data.frame(
    check = c("Missing", "Duplicated", "Unique IDs"),
    count = c(missing, duplicated_genes, unique_ids)
  ) %>%
    mutate(
      percentage = c(
        paste0(round((missing / total) * 100, 2), "%"),
        paste0(round((duplicated_genes / total) * 100, 2), "%"),
        paste0(round((unique_ids / total) * 100, 2), "%")
      )
    )
  dup_genes_from <- names(which(table(not_na1[[col_result_from]]) > 1))
  dup_genes_from <- dup_genes_from[dup_genes_from != ""]
  multiplicity_from <- round(sum(table(not_na1[[col_result_from]])[dup_genes_from])/
                               length(table(not_na1[[col_result_from]])[dup_genes_from]),2) 
  
  return(list(
    summary = summary,
    details = list(
      not_found = not_found,
      na_to = na_to,
      dup = duplicate,
      multiplicity = multiplicity_from)
  ))
}

# ============== MAPPATURA ATTRIBUTI ==============
get_biomart_attr <- function(id_type) {
  mapping <- list(
    SYMBOL = "external_gene_name",
    ENTREZID = "entrezgene_id",
    ENSEMBL = "ensembl_gene_id",
    UNIPROT = "uniprotswissprot",
    REFSEQ = "refseq_mrna",
    GENENAME = "description"
  )
  return(mapping[[id_type]])
}

get_gprofiler_target <- function(id_type) {
  mapping <- list(
    SYMBOL = "HGNC",
    ENTREZID = "ENTREZGENE_ACC",
    ENSEMBL = "ENSG",
    UNIPROT = "UNIPROTSWISSPROT"
  )
  return(mapping[[id_type]])
}

get_mygene_fields <- function(id_type) {
  mapping <- list(
    SYMBOL = "symbol",
    ENTREZID = "entrezgene",
    ENSEMBL = "ensembl.gene",
    UNIPROT = "uniprot.Swiss-Prot"
  )
  return(mapping[[id_type]])
}

extract_genes <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.data.frame(x)) return(x$gene)
  if (is.list(x) && "gene" %in% names(x)) return(x$gene)
  return(NA_character_)
}

# ============== FUNZIONI PER OGNI METODO ==============
run_biomart <- function(genes, species_dataset, from_type, to_type) {
  mart <- useEnsembl(biomart = "ensembl", dataset = species_dataset)
  attr_from <- get_biomart_attr(from_type)
  attr_to <- get_biomart_attr(to_type)

  if(from_type == "ENTREZID") {
    results <- getBM(attributes = c(attr_to, attr_from),
                     filters = attr_from,
                     values = genes,
                     mart = mart)
    colnames(results)[colnames(results) == attr_from] <- from_type
    colnames(results)[colnames(results) == attr_to] <- to_type
    results <- results[, c(from_type, to_type)]
  } else {
    results <- getBM(attributes = c(attr_from, attr_to),
                     filters = attr_from,
                     values = genes,
                     mart = mart)
    colnames(results)[colnames(results) == attr_from] <- from_type
    colnames(results)[colnames(results) == attr_to] <- to_type
    results <- results[, c(from_type, to_type)]
  }
  
  return(results)
}

run_mygene <- function(genes, species_code, from_type, to_type) {
  scopes <- tolower(get_mygene_fields(from_type))
  fields <- get_mygene_fields(to_type)
  
  if (species_code == "7227" && length(genes) > 300) {
    chunk_size <- 300
    chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
    
    gene_info_list <- vector("list", length(chunks))
    for (i in seq_along(chunks)) {
      cat(sprintf("mygene chunk %d/%d...\n", i, length(chunks)))
      tryCatch({
        gene_info_list[[i]] <- queryMany(chunks[[i]], scopes = scopes, 
                                         fields = fields, species = species_code)
        Sys.sleep(0.3)
      }, error = function(e) {
        cat(sprintf("  Chunk %d failed: %s\n", i, e$message))
        gene_info_list[[i]] <<- NULL
      })
    }
    gene_info_list <- Filter(Negate(is.null), gene_info_list)
    res1 <- data.frame()
    if(to_type == "ENSEMBL"){
      for (i in seq_along(gene_info_list)) {
        if(!is.null(gene_info_list[[i]]$ensembl.gene)){
      res <- data.frame(fromtp=gene_info_list[[i]]$query,
                        "ENSEMBL"=gene_info_list[[i]]$ensembl.gene)
      res1 <- rbind(res1,res)
        }
      }
      colnames(res1)[colnames(res1) == "fromtp"] <- from_type
      res <- res1[!is.na(res1[[to_type]]) & res1[[to_type]] != "",]
      res <- res[grep("^(ENS|WB|FB)", res[[to_type]]),]
    } else {
      for (i in seq_along(gene_info_list)) {
        if(!is.null(gene_info_list[[i]]$entrezgene) | !is.null(gene_info_list[[i]]$symbol)){
      res <- as.data.frame(gene_info_list[[i]]@listData)
      col_to <- switch(to_type,
                       ENTREZID = "X_id",
                       SYMBOL = "symbol",
                       "X_id"
      )
      res <- res[c("query", col_to)]
      res1 <- rbind(res1,res)
        }
      }
      res <- na.omit(res1)
      if(!col_to %in% colnames(res)) {
        warning(paste("Column", col_to, "not found in mygene results"))
        return(data.frame(setNames(list(character(0), character(0)), c(from_type, to_type))))
      }
    }
  } else {
    gene_info <- queryMany(genes, scopes = scopes, fields = fields, species = species_code)
  }
  
  if(species_code != "7227"){
    if(to_type == "ENSEMBL") {
      if(from_type == "ENTREZID" & (species_code == "10090" | species_code == "10116" | species_code == "6239")) {
        res <- data.frame("ENTREZID"=gene_info$query,
                          "ENSEMBL"=gene_info$ensembl.gene)
      res <- res[!is.na(res[[to_type]]) & res[[to_type]] != "",]
      res <- res[grep("^(ENS|WB|FB)", res[[to_type]]),]
      } else if(from_type == "ENTREZID" & (species_code == "7955" | species_code == "human")) {
        queries <- gene_info$query
        ensembls <- gene_info$ensembl
        res <- do.call(rbind, Map(function(q, e) {
          genes <- extract_genes(e)
          data.frame(ENTREZID = q, ENSEMBL = genes)
        }, queries, ensembls))
        rownames(res) <- NULL
          res <- res[!is.na(res[[to_type]]) & res[[to_type]] != "",]
          res <- res[grep("^(ENS|WB|FB)", res[[to_type]]),]
      } else if(from_type == "SYMBOL" &  species_code == "human") {
        queries <- gene_info$query
        ensembls <- gene_info$ensembl
        res <- do.call(rbind, Map(function(q, e) {
          genes <- extract_genes(e)
          data.frame(ENTREZID = q, ENSEMBL = genes)
        }, queries, ensembls))
        rownames(res) <- NULL
        res <- res[!is.na(res[[to_type]]) & res[[to_type]] != "",]
        res <- res[grep("^(ENS|WB|FB)", res[[to_type]]),]
      } else if(from_type == "SYMBOL" & (species_code == "3702" | species_code == "10116" | species_code == "6239")){ 
        res <- data.frame("SYMBOL"=gene_info$query,
                          "ENSEMBL"=gene_info$ensembl.gene)
        res <- res[!is.na(res[[to_type]]) & res[[to_type]] != "",]
        res <- res[grep("^(ENS|WB|FB)", res[[to_type]]),]
      } else {
        res <- data.frame(
        FROM = gene_info@listData$query,
        TO = gene_info@listData$`_id`
      )
      colnames(res) <- c(from_type, to_type)
      res <- res[!is.na(res[[to_type]]) & res[[to_type]] != "",]
      res <- res[grep("^(ENS|WB|FB)", res[[to_type]]),]
      }} else {
      res <- as.data.frame(gene_info@listData)
    
    col_to <- switch(to_type,
                     ENTREZID = "X_id",
                     SYMBOL = "symbol",
                     "X_id"
    )
    res <- res[c("query", col_to)]
    res <- na.omit(res)
    if(!col_to %in% colnames(res)) {
      warning(paste("Column", col_to, "not found in mygene results"))
      return(data.frame(setNames(list(character(0), character(0)), c(from_type, to_type))))
    }
      }
  }
  colnames(res) <- c(from_type, to_type)
  if(to_type == "ENTREZID" | to_type == "SYMBOL") {
    res <- res[!grepl("^(ENS|WB|FB)", res[[to_type]]),]
  }
  
  return(res)
}

run_orgdb <- function(genes, orgdb, from_type, to_type) {
  available_keys <- keytypes(orgdb)
  from_key <- from_type
  to_key <- to_type
  
  if(from_type == "ENSEMBL" && !"ENSEMBL" %in% available_keys && "TAIR" %in% available_keys) {
    from_key <- "TAIR"
  }
  if(to_type == "ENSEMBL" && !"ENSEMBL" %in% available_keys && "TAIR" %in% available_keys) {
    to_key <- "TAIR"
  }
  
  if(from_type == "SYMBOL" && !"SYMBOL" %in% available_keys && "GENENAME" %in% available_keys) {
    from_key <- "GENENAME"
  }
  if(to_type == "SYMBOL" && !"SYMBOL" %in% available_keys && "GENENAME" %in% available_keys) {
    to_key <- "GENENAME"
  }
  
  res <- AnnotationDbi::select(orgdb, keys = as.character(genes), columns = to_key, keytype = from_key)
  colnames(res) <- c(from_type, to_type)
  return(res)
}

run_gprofiler <- function(genes, organism, from_type, to_type) {
 
   if (to_type == "SYMBOL"){
    if (organism == "hsapiens"){
    target <- get_gprofiler_target(to_type)
  } else if (organism == "mmusculus"){
    target <- "MGI"
  } else if (organism == "drerio"){
    target <- "ZFIN_ID"
  } else if (organism == "scerevisiae"){
    target <- "ENSG"
  } else if (organism == "athaliana"){
    target <- "TAIR_LOCUS"
  }else if (organism == "rnorvegicus"){
    target <- "RGD"
  }else if (organism == "dmelanogaster"){
    target <- "FLYBASENAME_GENE"
  }else if (organism == "celegans"){
    target <- "WORMBASE_LOCUS"
  }} else {
    target <- get_gprofiler_target(to_type)
  }
  
  
  if(from_type == "ENTREZID") {
    if(organism == "dmelanogaster"){
      res <- gconvert(as.character(genes), organism = organism, 
                      target = target, 
                      numeric_ns = "ENTREZGENE",
                      mthreshold = Inf, filter_na = TRUE)
      res <- res[, c("input", "target")]
    } else {
    res <- gconvert(as.character(genes), organism = organism, target = target, 
                    numeric_ns = "ENTREZGENE_ACC",
                    mthreshold = Inf, filter_na = TRUE)
    res <- res[, c("input", "target")]
    }
  } else if(to_type == "ENTREZID" & organism == "dmelanogaster") {
    res <- gconvert(as.character(genes), organism = organism, 
                    target = "ENTREZGENE", 
                    mthreshold = Inf, filter_na = TRUE)
    res <- res[, c("input", "target")]
  } else {
  res <- gconvert(as.character(genes), organism = organism, target = target, 
                  mthreshold = Inf, filter_na = TRUE)
  res <- res[, c("input", "target")]
  } 
  colnames(res) <- c(from_type, to_type)
  res <- res %>%
    mutate(!!to_type := ifelse(grepl("^None", .data[[to_type]]), NA, .data[[to_type]]))
  
  return(res)
}


run_geneslator <- function(genes, orgdb, from_type, to_type) {
  res <- geneslator::select(orgdb, keys = as.character(genes), 
                            columns = c(from_type, to_type), keytype = from_type)
  return(res)
}

# ============== JACCARD + FISHER TEST ==============
calc_jaccard_fisher <- function(df1, df2, col, universo = NULL) {
  
  vals1 <- df1[[col]][!is.na(df1[[col]]) & df1[[col]] != ""]
  vals2 <- df2[[col]][!is.na(df2[[col]]) & df2[[col]] != ""]
  
  a <- length(intersect(vals1, vals2))
  b <- length(setdiff(vals1, vals2))
  c <- length(setdiff(vals2, vals1))
  
  if(!is.null(universo)) {
    d <- length(setdiff(universo, union(vals1, vals2)))
  } else {
    d <- 0
  }
  
  jaccard_idx <- ifelse((a + b + c) > 0, a / (a + b + c), NA)
  
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  fisher_res <- fisher.test(tab)
  
  return(list(
    jaccard_index = jaccard_idx,
    pvalue = fisher_res$p.value,
    odds_ratio = as.numeric(fisher_res$estimate),
    contingency = list(a = a, b = b, c = c, d = d)
  ))
}

# ============== PIPELINE ==============
run_comparison_pipeline <- function(input_genes,
                                    from_type,
                                    to_type,
                                    species_dataset,
                                    species_code,
                                    organism_gprofiler,
                                    orgdb,
                                    orgdb_geneslator,
                                    reference = "Geneslator",
                                    universo = NULL,
                                    output_file = "comparison_results.xlsx",
                                    methods_to_run = c("Geneslator", "Biomart", "Org.Db", "mygene", "gprofiler2")) {
  
  if(is.data.frame(input_genes)) {
    genes <- as.character(input_genes[[1]])
  } else {
    genes <- as.character(input_genes)
  }
  
  genes <- genes[!is.na(genes) & genes != ""]
  total <- length(genes)
  
  message(paste("Converting", total, "genes from", from_type, "to", to_type))
  
  method_functions <- list(
    Geneslator = function() run_geneslator(genes, orgdb_geneslator, from_type, to_type),
    Biomart = function() run_biomart(genes, species_dataset, from_type, to_type),
    Org.Db = function() run_orgdb(genes, orgdb, from_type, to_type),
    mygene = function() run_mygene(genes, species_code, from_type, to_type),
    gprofiler2 = function() run_gprofiler(genes, organism_gprofiler, from_type, to_type)
  )
  
  methods_results <- list()
  for(m in methods_to_run) {
    message(paste("Running", m, "..."))
    tryCatch({
      methods_results[[m]] <- method_functions[[m]]()
    }, error = function(e) {
      warning(paste("Error in", m, ":", e$message))
      methods_results[[m]] <<- data.frame(setNames(list(character(0), character(0)), c(from_type, to_type)))
    })
  }
  
  message("Checking conversions...")

  checks <- lapply(methods_results, function(res) {
    check_conversion(
      example = data.frame(setNames(list(genes), from_type)),
      result = res,
      col_example = from_type,
      col_result_from = from_type,
      col_result_to = to_type
    )
  })
  
  message("Building summary table...")
  
  get_metric <- function(check_res, metric) {
    idx <- which(check_res$summary$check == metric)
    return(check_res$summary$percentage[idx])
  }
  
  get_multiplicity <- function(check_res) {
    mult <- check_res$details$multiplicity
    if(is.na(mult) || is.null(mult)) return("-")
    return(paste0("(", mult, ")"))
  }

  message("Calculating Jaccard index and Fisher p-values...")
  ref_df <- methods_results[[reference]]
  method_names <- names(methods_results)
  
  jaccard_results <- lapply(method_names, function(m) {
    if(m == reference) {
      return(list(jaccard_index = NA, pvalue = NA))
    }
    tryCatch({
      calc_jaccard_fisher(ref_df, methods_results[[m]], col = to_type, universo = universo)
    }, error = function(e) {
      warning(paste("Error comparing", reference, "vs", m, ":", e$message))
      return(list(jaccard_index = NA, pvalue = NA))
    })
  })
  names(jaccard_results) <- method_names
  
  p_raw <- sapply(method_names, function(m) {
    if(m == reference) return(NA)
    pval <- jaccard_results[[m]]$pvalue
    if(is.na(pval)) return(NA)
    return(pval)
  })
  
  p_adjusted <- rep(NA, length(p_raw))
  valid <- !is.na(p_raw) & method_names != reference
  p_adjusted[valid] <- p.adjust(p_raw[valid], method = "fdr")
  
  pvalues <- sapply(seq_along(method_names), function(i) {
    m <- method_names[i]
    if(m == reference) return("-")
    pval <- p_adjusted[i]
    if(is.na(pval)) return("NA")
    if(pval < 0.001) return("<0.001")
    return(round(pval, 3))
  })
  names(pvalues) <- method_names
  
  jaccard_idx <- sapply(method_names, function(m) {
    if(m == reference) return("-")
    jidx <- jaccard_results[[m]]$jaccard_index
    if(is.na(jidx)) return("NA")
    return(round(jidx, 3))
  })
  
  summary_table <- data.frame(
    Conversion = c("One to one", "Unmapped", "One to many", "(Multiplicity)", 
                   "Jaccard Index", "Fisher p-value")
  )
  
  for(m in method_names) {
    check_res <- checks[[m]]
    summary_table[[m]] <- c(
      get_metric(check_res, "Unique IDs"),
      get_metric(check_res, "Missing"),
      get_metric(check_res, "Duplicated"),
      get_multiplicity(check_res),
      jaccard_idx[m],
      pvalues[m]
    )
  }

  message("Saving to Excel...")
  write.xlsx(summary_table, output_file, rowNames = FALSE)
  message(paste("Done! Results saved to:", output_file))
  
  return(list(
    summary_table = summary_table,
    methods_results = methods_results,
    checks = checks,
    jaccard_results = jaccard_results
  ))
}

setwd("./")

#Homo sapiens
GeneslatorDb("Homo sapiens")
#SYMBOL
human <- readRDS("Human/example.symbolid.RDS") 
human <- unique(human)
colnames(human) <- "SYMBOL"

results_hs <- run_comparison_pipeline(
  input_genes = human$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "hsapiens_gene_ensembl",
  species_code = "human",
  organism_gprofiler = "hsapiens",
  orgdb = org.Hs.eg.db,
  orgdb_geneslator = org.Hsapiens.db,
  output_file = "human_symbol_to_ensembl_march26.xlsx"
)

results_hs <- run_comparison_pipeline(
  input_genes = human$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "hsapiens_gene_ensembl",
  species_code = "human",
  organism_gprofiler = "hsapiens",
  orgdb = org.Hs.eg.db,
  orgdb_geneslator = org.Hsapiens.db,
  output_file = "human_symbol_to_entrez_march26.xlsx"
)
#ENTREZ
human <- readRDS("example.entrezid_human.RDS") #pulire gli AC, NC, RA
human <- unique(human)
colnames(human) <- "ENTREZID"

results_hs <- run_comparison_pipeline(
  input_genes = human$ENTREZID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "hsapiens_gene_ensembl",
  species_code = "human",
  organism_gprofiler = "hsapiens",
  orgdb = org.Hs.eg.db,
  orgdb_geneslator = org.Hsapiens.db,
  output_file = "human_entrez_to_ensembl_march26.xlsx"
)

results_hs <- run_comparison_pipeline(
  input_genes = human$ENTREZID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "hsapiens_gene_ensembl",
  species_code = "human",
  organism_gprofiler = "hsapiens",
  orgdb = org.Hs.eg.db,
  orgdb_geneslator = org.Hsapiens.db,
  output_file = "human_entrez_to_symbol_march26.xlsx"
)

#ENSEMBL
human <- readRDS("list.ensembl.RDS") 
human <- unique(human)
colnames(human) <- "ENSEMBL"


results_hs <- run_comparison_pipeline(
  input_genes = human$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "hsapiens_gene_ensembl",
  species_code = "human",
  organism_gprofiler = "hsapiens",
  orgdb = org.Hs.eg.db,
  orgdb_geneslator = org.Hsapiens.db,
  output_file = "human_ensembl_to_entrezid_march26.xlsx"
)

results_hs <- run_comparison_pipeline(
  input_genes = human$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "hsapiens_gene_ensembl",
  species_code = "human",
  organism_gprofiler = "hsapiens",
  orgdb = org.Hs.eg.db,
  orgdb_geneslator = org.Hsapiens.db,
  output_file = "human_ensembl_to_symbol_march26.xlsx"
)

#Mouse 
GeneslatorDb("Mus musculus")
#SYMBOL
mouse1 <- readRDS("Mouse/example.symbolid.RDS") #pulire gli AC, NC, RA
mouse1 <- unique(mouse1)
rownames(mouse1) <- NULL
mouse1 <- as.data.frame(mouse1[-29759,])
rownames(mouse1) <- NULL
colnames(mouse1) <- "SYMBOL"


results_mm <- run_comparison_pipeline(
  input_genes = mouse1$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "mmusculus_gene_ensembl",
  species_code = "10090",
  organism_gprofiler = "mmusculus",
  orgdb = org.Mm.eg.db,
  orgdb_geneslator = org.Mmusculus.db,
  output_file = "Mouse_symbol_to_ensembl_march26.xlsx"
)

results_mm <- run_comparison_pipeline(
  input_genes = mouse1$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "mmusculus_gene_ensembl",
  species_code = "10090",
  organism_gprofiler = "mmusculus",
  orgdb = org.Mm.eg.db,
  orgdb_geneslator = org.Mmusculus.db,
  output_file = "Mouse_symbol_to_entrezid_march26.xlsx"
)

#ENSEMBL
mouse1 <- readRDS("Mouse/example.ensembl.RDS")
mouse1$Ensembl <- as.character(mouse1$Ensembl)

results_mm <- run_comparison_pipeline(
  input_genes = mouse1$Ensembl,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "mmusculus_gene_ensembl",
  species_code = "10090",
  organism_gprofiler = "mmusculus",
  orgdb = org.Mm.eg.db,
  orgdb_geneslator = org.Mmusculus.db,
  output_file = "Mouse_ensembl_to_symbol_march26.xlsx"
)

results_mm <- run_comparison_pipeline(
  input_genes = mouse1$Ensembl,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "mmusculus_gene_ensembl",
  species_code = "10090",
  organism_gprofiler = "mmusculus",
  orgdb = org.Mm.eg.db,
  orgdb_geneslator = org.Mmusculus.db,
  output_file = "Mouse_ensembl_to_entrez_march26.xlsx"
)

#ENTREZID
mouse1 <- readRDS("Mouse/example.entrezid.RDS")
mouse1 <- unique(mouse1)
rownames(mouse1) <- NULL
mouse1$GeneID <- as.character(mouse1$GeneID)

results_mm <- run_comparison_pipeline(
  input_genes = mouse1$GeneID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "mmusculus_gene_ensembl",
  species_code = "10090",
  organism_gprofiler = "mmusculus",
  orgdb = org.Mm.eg.db,
  orgdb_geneslator = org.Mmusculus.db,
  output_file = "Mouse_entrez_to_symbol_march26.xlsx"
)

results_mm <- run_comparison_pipeline(
  input_genes = mouse1$GeneID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "mmusculus_gene_ensembl",
  species_code = "10090",
  organism_gprofiler = "mmusculus",
  orgdb = org.Mm.eg.db,
  orgdb_geneslator = org.Mmusculus.db,
  output_file = "Mouse_entrezid_to_ensembl_march26.xlsx"
)


#Zebrafish
GeneslatorDb("Danio rerio")
#SYMBOL
input_df <- read.csv("Zebrafish/GSE306907_gene_count_matrix.csv/conte_zebrafish_symbol.txt", sep = ",", header = FALSE)
colnames(input_df)[1] <- "SYMBOL"
input_df <- as.data.frame(input_df$SYMBOL)
colnames(input_df)[1] <- "SYMBOL"
input_df$SYMBOL <- gsub("_",":",input_df$SYMBOL)


results_zf <- run_comparison_pipeline(
  input_genes = input_df$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "drerio_gene_ensembl",
  species_code = "7955",
  organism_gprofiler = "drerio",
  orgdb = org.Dr.eg.db,
  orgdb_geneslator = org.Drerio.db,
  output_file = "zebrafish_symbol_to_ensembl_march26.xlsx"
)


results_zf <- run_comparison_pipeline(
  input_genes = input_df$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "drerio_gene_ensembl",
  species_code = "7955",
  organism_gprofiler = "drerio",
  orgdb = org.Dr.eg.db,
  orgdb_geneslator = org.Drerio.db,
  output_file = "zebrafish_symbol_to_entrezid_march26.xlsx"
)

#ENSEMBL
zebr <- read_excel("Zebrafish/GSE274820_Repository_AnalysisData_2019_A5806_DESeq.xlsx", sheet = 3) 
colnames(zebr)[1] <- "ENSEMBL"
zebr <- as.data.frame(zebr$ENSEMBL)
colnames(zebr)[1] <- "ENSEMBL"

results_zf <- run_comparison_pipeline(
  input_genes = zebr$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "drerio_gene_ensembl",
  species_code = "7955",
  organism_gprofiler = "drerio",
  orgdb = org.Dr.eg.db,
  orgdb_geneslator = org.Drerio.db,
  output_file = "zebrafish_ensembl_to_symbol_march26.xlsx"
)

results_zf <- run_comparison_pipeline(
  input_genes = zebr$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "drerio_gene_ensembl",
  species_code = "7955",
  organism_gprofiler = "drerio",
  orgdb = org.Dr.eg.db,
  orgdb_geneslator = org.Drerio.db,
  output_file = "zebrafish_ensembl_to_entrezid_march26.xlsx"
)

#ENTREZID
zeb1 <- read.csv("Zebrafish/Danio_rerio_EntrezID_GPL18967.txt", sep="\t") 
colnames(zeb1) <- "gene_name"
zeb1$gene_name <- as.character(zeb1$gene_name)

results_zf <- run_comparison_pipeline(
  input_genes = zeb1$gene_name,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "drerio_gene_ensembl",
  species_code = "7955",
  organism_gprofiler = "drerio",
  orgdb = org.Dr.eg.db,
  orgdb_geneslator = org.Drerio.db,
  output_file = "zebrafish_entrezid_to_symbol_march26.xlsx"
)

results_zf <- run_comparison_pipeline(
  input_genes = zeb1$gene_name,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "drerio_gene_ensembl",
  species_code = "7955",
  organism_gprofiler = "drerio",
  orgdb = org.Dr.eg.db,
  orgdb_geneslator = org.Drerio.db,
  output_file = "zebrafish_entrezid_to_ensembl_march26.xlsx"
)


#Yeast
GeneslatorDb("Saccharomyces cerevisiae")
#SYMBOL
yeast1 <- read.csv("Yeast/GSE280426_Normalize_Counts.txt/lista_geni_yeast.txt", 
                   sep = ";", header=FALSE)
yeast1$V1 <- gsub("_", "", yeast1$V1)
yeast1 <- as.data.frame(yeast1[,c(1)])
yeast1 <- unique(yeast1)
colnames(yeast1) <- "SYMBOL"

results <- run_comparison_pipeline(
  input_genes = yeast1$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "scerevisiae_gene_ensembl",
  species_code = "4932",
  organism_gprofiler = "scerevisiae",
  orgdb = org.Sc.sgd.db,
  orgdb_geneslator = org.Scerevisiae.db,
  output_file = "scerevisiae_symbol_entrezid_march26.xlsx"
)

results <- run_comparison_pipeline(
  input_genes = yeast1$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "scerevisiae_gene_ensembl",
  species_code = "4932",
  organism_gprofiler = "scerevisiae",
  orgdb = org.Sc.sgd.db,
  orgdb_geneslator = org.Scerevisiae.db,
  output_file = "scerevisiae_symbol_to_ensembl_march26.xlsx"
)

#ENTREZID 
scer1 <- read.csv("Yeast/Saccharomyces_cerevisiae_ENTREZID_GPL90.txt", sep = "\t", header = FALSE) 
colnames(scer1) <- "GeneID"
scer1$GeneID <- as.character(scer1$GeneID)

results <- run_comparison_pipeline(
  input_genes = scer1$GeneID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "scerevisiae_gene_ensembl",
  species_code = "4932",
  organism_gprofiler = "scerevisiae",
  orgdb = org.Sc.sgd.db,
  orgdb_geneslator = org.Scerevisiae.db,
  output_file = "scerevisiae_entrezid_to_symbol_march26.xlsx"
)

results <- run_comparison_pipeline(
  input_genes = scer1$GeneID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "scerevisiae_gene_ensembl",
  species_code = "4932",
  organism_gprofiler = "scerevisiae",
  orgdb = org.Sc.sgd.db,
  orgdb_geneslator = org.Scerevisiae.db,
  output_file = "scerevisiae_entrezid_to_ensembl_march26.xlsx"
)

#ENSEMBL 
yeast1 <- read.csv("Yeast/GSE307461_raw_gene_counts.tsv/GSE307461_raw_gene_counts.tsv", sep="\t") #pulire gli AC, NC, RA
yeast1 <- as.data.frame(yeast1$Geneid)
colnames(yeast1) <- "Geneid"

results <- run_comparison_pipeline(
  input_genes = yeast1$Geneid,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "scerevisiae_gene_ensembl",
  species_code = "4932",
  organism_gprofiler = "scerevisiae",
  orgdb = org.Sc.sgd.db,
  orgdb_geneslator = org.Scerevisiae.db,
  output_file = "scerevisiae_ensembl_to_symbol_march26.xlsx"
)

results <- run_comparison_pipeline(
  input_genes = yeast1$Geneid,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "scerevisiae_gene_ensembl",
  species_code = "4932",
  organism_gprofiler = "scerevisiae",
  orgdb = org.Sc.sgd.db,
  orgdb_geneslator = org.Scerevisiae.db,
  output_file = "scerevisiae_ensembl_to_entrezid_march26.xlsx"
)

#Arabidopsis
GeneslatorDb("Arabidopsis thaliana")
#SYMBOL
ara1 <- read.csv("Arabidopsis/GSM8858406_TPM_Ler_R1.csv", sep=",") 
ara1 <- as.data.frame(ara1$Ara_gene)
ara1$`ara1$Ara_gene` <- gsub("\\|.*", "", ara1$`ara1$Ara_gene`)
colnames(ara1) <- "Gene"
ara1 <- unique(ara1)

results <- run_comparison_pipeline(
  input_genes = ara1$Gene,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "athaliana_gene_ensembl",
  species_code = "3702",
  organism_gprofiler = "athaliana",
  orgdb = org.At.tair.db,
  orgdb_geneslator = org.Athaliana.db,
  output_file = "arabidopsis_symbol_entrezid_march26.xlsx"
)

results <- run_comparison_pipeline(
  input_genes = ara1$Gene,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "athaliana_gene_ensembl",
  species_code = "3702",
  organism_gprofiler = "athaliana",
  orgdb = org.At.tair.db,
  orgdb_geneslator = org.Athaliana.db,
  output_file = "arabidopsis_symbol_to_ensembl_march26.xlsx"
)

#ENSEMBL
arab <- read.csv("Arabidopsis/GSE292906_RNAseq-count_all.txt", sep="\t") 
arab <- as.data.frame(arab$Geneid)
colnames(arab) <- "Geneid"

results <- run_comparison_pipeline(
  input_genes = arab$Geneid,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "athaliana_gene_ensembl",
  species_code = "3702",
  organism_gprofiler = "athaliana",
  orgdb = org.At.tair.db,
  orgdb_geneslator = org.Athaliana.db,
  output_file = "arabidopsis_ensembl_to_symbol_march26.xlsx"
)

results <- run_comparison_pipeline(
  input_genes = arab$Geneid,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "athaliana_gene_ensembl",
  species_code = "3702",
  organism_gprofiler = "athaliana",
  orgdb = org.At.tair.db,
  orgdb_geneslator = org.Athaliana.db,
  output_file = "arabidopsis_ensembl_to_entrez_march26.xlsx"
)

#ENTREZID
ara1 <- read.csv("Arabidopsis/Arabidopsis_thaliana_EntrezID_GPL198.txt", sep="\t", header= FALSE) 
colnames(ara1) <- "GeneID"
ara1$GeneID <- as.character(ara1$GeneID)

results <- run_comparison_pipeline(
  input_genes = ara1$GeneID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "athaliana_gene_ensembl",
  species_code = "3702",
  organism_gprofiler = "athaliana",
  orgdb = org.At.tair.db,
  orgdb_geneslator = org.Athaliana.db,
  output_file = "arabidopsis_entrez_to_symbol_march26.xlsx"
)

results <- run_comparison_pipeline(
  input_genes = ara1$GeneID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "athaliana_gene_ensembl",
  species_code = "3702",
  organism_gprofiler = "athaliana",
  orgdb = org.At.tair.db,
  orgdb_geneslator = org.Athaliana.db,
  output_file = "arabidopsis_entrez_to_ensembl_march26.xlsx"
)


#Rattus norvegicus
GeneslatorDb("Rattus norvegicus")
#SYMBOL
rattus <- readRDS("Rattus/example.symbolid.RDS") 
rattus <- unique(rattus)
colnames(rattus) <- "SYMBOL"

results_rn <- run_comparison_pipeline(
  input_genes = rattus$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "rnorvegicus_gene_ensembl",
  species_code = "10116",
  organism_gprofiler = "rnorvegicus",
  orgdb = org.Rn.eg.db,
  orgdb_geneslator = org.Rnorvegicus.db,
  output_file = "rattus_symbol_to_ensembl_march26.xlsx"
)

results_rn <- run_comparison_pipeline(
  input_genes = rattus$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "rnorvegicus_gene_ensembl",
  species_code = "10116",
  organism_gprofiler = "rnorvegicus",
  orgdb = org.Rn.eg.db,
  orgdb_geneslator = org.Rnorvegicus.db,
  output_file = "rattus_symbol_to_entrez_march26.xlsx"
)

#ENTREZ
rattus <- readRDS("Rattus/entrezid.RDS") 
rattus <- unique(rattus)
colnames(rattus) <- "ENTREZID"

results_rn <- run_comparison_pipeline(
  input_genes = rattus$ENTREZID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "rnorvegicus_gene_ensembl",
  species_code = "10116",
  organism_gprofiler = "rnorvegicus",
  orgdb = org.Rn.eg.db,
  orgdb_geneslator = org.Rnorvegicus.db,
  output_file = "rattus_entrez_to_ensembl_march26.xlsx"
)


results_rn <- run_comparison_pipeline(
  input_genes = rattus$ENTREZID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "rnorvegicus_gene_ensembl",
  species_code = "10116",
  organism_gprofiler = "rnorvegicus",
  orgdb = org.Rn.eg.db,
  orgdb_geneslator = org.Rnorvegicus.db,
  output_file = "rattus_entrez_to_symbol_march26.xlsx"
)

#ENSEMBL
rattus <- readRDS("Rattus/ensembl.RDS") #pulire gli AC, NC, RA
rattus <- unique(rattus)
colnames(rattus) <- "ENSEMBL"

results_rn <- run_comparison_pipeline(
  input_genes = rattus$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "rnorvegicus_gene_ensembl",
  species_code = "10116",
  organism_gprofiler = "rnorvegicus",
  orgdb = org.Rn.eg.db,
  orgdb_geneslator = org.Rnorvegicus.db,
  output_file = "rattus_ensembl_to_entrezid_march26.xlsx"
)

results_rn <- run_comparison_pipeline(
  input_genes = rattus$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "rnorvegicus_gene_ensembl",
  species_code = "10116",
  organism_gprofiler = "rnorvegicus",
  orgdb = org.Rn.eg.db,
  orgdb_geneslator = org.Rnorvegicus.db,
  output_file = "rattus_ensembl_to_symbol_march26.xlsx"
)

#Drosophila
GeneslatorDb("Drosophila melanogaster")
#SYMBOL
input_df <- readRDS("Drosophila/example.symbolid.RDS")
colnames(input_df)[1] <- "SYMBOL"
input_df <- unique(input_df)

results_dm_symbol_entrezid <- run_comparison_pipeline(
  input_genes = input_df$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "dmelanogaster_gene_ensembl",
  species_code = "7227",
  organism_gprofiler = "dmelanogaster",
  orgdb = org.Dm.eg.db,
  orgdb_geneslator = org.Dmelanogaster.db,
  output_file = "drosophila_symbol_entrezid_march26.xlsx"
)

results_dm_symbol_ensembl <- run_comparison_pipeline(
  input_genes = input_df$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "dmelanogaster_gene_ensembl",
  species_code = "7227",
  organism_gprofiler = "dmelanogaster",
  orgdb = org.Dm.eg.db,
  orgdb_geneslator = org.Dmelanogaster.db,
  output_file = "drosophila_symbol_ensembl_march26.xlsx"
)

#ENSEMBL
input_df <- readRDS("Drosophila/ensembl.RDS")
colnames(input_df)[1] <- "ENSEMBL"
input_df<- unique(input_df)

results_dm_ensembl_entrezid <- run_comparison_pipeline(
  input_genes = input_df$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "dmelanogaster_gene_ensembl",
  species_code = "7227",
  organism_gprofiler = "dmelanogaster",
  orgdb = org.Dm.eg.db,
  orgdb_geneslator = org.Dmelanogaster.db,
  output_file = "drosophila_ensembl_entrezid_march26.xlsx"
)

results_dm_ensembl_symbol <- run_comparison_pipeline(
  input_genes = input_df$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "dmelanogaster_gene_ensembl",
  species_code = "7227",
  organism_gprofiler = "dmelanogaster",
  orgdb = org.Dm.eg.db,
  orgdb_geneslator = org.Dmelanogaster.db,
  output_file = "drosophila_esembl_symbol_march26.xlsx"
)


#ENTREZID
input_df <- readRDS("Drosophila/entrez.RDS")
colnames(input_df)[1] <- "ENTREZID"
input_df <- unique(input_df)

results_dm_entrezid_ensembl<- run_comparison_pipeline(
  input_gene = input_df$ENTREZID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "dmelanogaster_gene_ensembl",
  species_code = "7227",
  organism_gprofiler = "dmelanogaster",
  orgdb = org.Dm.eg.db,
  orgdb_geneslator = org.Dmelanogaster.db,
  output_file = "drosophila_entrezid_ensembl_march26.xlsx"
)

results_dm_entrezid_symbol <- run_comparison_pipeline(
  input_gene = input_df$ENTREZID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "dmelanogaster_gene_ensembl",
  species_code = "7227",
  organism_gprofiler = "dmelanogaster",
  orgdb = org.Dm.eg.db,
  orgdb_geneslator = org.Dmelanogaster.db,
  output_file = "drosophila_entrez_symbol_march26.xlsx"
)

#Celegans
GeneslatorDb("Caenorhabditis elegans")
#SYMBOL
input_df <- readRDS("Celegans/example.symbolid.RDS")
colnames(input_df)[1] <- "SYMBOL"
input_df <- unique(input_df)

results_ce_symbol_entrezid <- run_comparison_pipeline(
  input_genes = input_df$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENTREZID",
  species_dataset = "celegans_gene_ensembl",
  species_code = "6239",
  organism_gprofiler = "celegans",
  orgdb = org.Ce.eg.db,
  orgdb_geneslator = org.Celegans.db,
  output_file = "celegans_symbol_entrezid_march26.xlsx" #mygene da sistemare
)

results_ce_symbol_ensembl <- run_comparison_pipeline(
  input_gene = input_df$SYMBOL,
  from_type = "SYMBOL",
  to_type = "ENSEMBL",
  species_dataset = "celegans_gene_ensembl",
  species_code = "6239",
  organism_gprofiler = "celegans",
  orgdb = org.Ce.eg.db,
  orgdb_geneslator = org.Celegans.db,
  output_file = "celegans_symbol_ensembl_march26.xlsx"
)


#ENSEMBL
input_df <- readRDS("Celegans/ensembl.RDS")
colnames(input_df)[1] <- "ENSEMBL"
input_df <- unique(input_df)

results_ce_ensembl_entrezid <- run_comparison_pipeline(
  input_genes = input_df$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  species_dataset = "celegans_gene_ensembl",
  species_code = "6239",
  organism_gprofiler = "celegans",
  orgdb = org.Ce.eg.db,
  orgdb_geneslator = org.Celegans.db,
  output_file = "celegans_symbol_entrezid_march26.xlsx" 
)

results_ce_ensembl_symbol <- run_comparison_pipeline(
  input_genes = input_df$ENSEMBL,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  species_dataset = "celegans_gene_ensembl",
  species_code = "6239",
  organism_gprofiler = "celegans",
  orgdb = org.Ce.eg.db,
  orgdb_geneslator = org.Celegans.db,
  output_file = "celegans_ensembl_symbol_march26.xlsx"
)


#ENTREZID
input_df <- readRDS("Celegans/entrezid.RDS")
colnames(input_df)[1] <- "ENTREZID"
input_df <- unique(input_df)

results_ce_entrezid_ensembl <- run_comparison_pipeline(
  input_genes = input_df$ENTREZID,
  from_type = "ENTREZID",
  to_type = "ENSEMBL",
  species_dataset = "celegans_gene_ensembl",
  species_code = "6239",
  organism_gprofiler = "celegans",
  orgdb = org.Ce.eg.db,
  orgdb_geneslator = org.Celegans.db,
  output_file = "celegans_entrezid_ensembl_march26.xlsx" #mygene da sistemare
)

results_ce_entrez_symbol <- run_comparison_pipeline(
  input_genes = input_df$ENTREZID,
  from_type = "ENTREZID",
  to_type = "SYMBOL",
  species_dataset = "celegans_gene_ensembl",
  species_code = "6239",
  organism_gprofiler = "celegans",
  orgdb = org.Ce.eg.db,
  orgdb_geneslator = org.Celegans.db,
  output_file = "celegans_entrez_symbol_march26.xlsx"
)


