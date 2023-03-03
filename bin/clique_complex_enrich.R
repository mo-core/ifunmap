#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
library(argparse)
library(tools)
library(igraph)
suppressPackageStartupMessages(library(WebGestaltR))

# for each clique, we do ORA analysis with respect to
# CORUM, BioPlex, GO_BP, GO_MF, GO_CC
# for each database, we obtain the top enriched item
# bases on its FDR, we classify the clique into 4 categories
# (1) CORUM FDR < 0.05  (class 1)
# (2) CORUM >= 0.05, BioPlex < 0.05 (class 2)
# (3) COURM >= 0.05, BioPlex >= 0.05, GO_BP or GO_MF < 0.05 and GO_CC < 0.05 (class 3)
# (4) others (class 4)

read_gmt <- function(file){
  if(!grepl("\\.gmt$",file)[1]){stop("Pathway information must be a .gmt file")}
  geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
  geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
  names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
  geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
  geneSetDB = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
  return(geneSetDB)
}

parser <- ArgumentParser()
parser$add_argument("--clique-file", type = "character",
    required = TRUE,
    help = "path to the output file from ICE")
parser$add_argument("--mapping-file", type = "character",
    required = TRUE,
    help = "path to the gene id mapping file")
parser$add_argument("--corum-gmt", type = "character",
    required = TRUE,
    help = "path to the corum gmt file")
parser$add_argument("--bioplex-gmt", type = "character",
    required = TRUE,
    help = "path to the bioplex gmt file")
parser$add_argument("--network-file", type = "character",
    required = TRUE,
    help = "path to the network edge list file")

args <- parser$parse_args()

input_clique_file <- args$clique_file
network_edge_list_file <- args$network_file
corum_gmt_file <- args$corum_gmt
bioplex_gmt_file <- args$bioplex_gmt
id_mapping_file <- args$mapping_file
sigmethod <- "top"
cutoff <- 1

# create gene symbol to ensembl_gene_id mapping table
id_mapping <- read.table(id_mapping_file, header=TRUE, sep="\t")
symbol_to_id <- id_mapping$ensembl_gene_id
names(symbol_to_id) <- id_mapping$gene_symbol

# create a union graph of BioPlex, BioGrid, HI-Union
bioplex_edgelist_file <- "./BioPlex_3_edge_list.tsv"
bioplex_network <- read_graph(bioplex_edgelist_file, format="ncol", directed=FALSE)

biogrid_edglist_file <- "./BioGRID_BioGRID.tsv"
biogrid_network <- read_graph(biogrid_edglist_file, format="ncol", directed=FALSE)

hi_uion_edgelist_file <- "./HI-union.tsv"
hi_uion_network <- read_graph(hi_uion_edgelist_file, format="ncol", directed=FALSE)

combined_network <- bioplex_network %u% biogrid_network %u% hi_uion_network
combined_network <- simplify(combined_network)
stopifnot(is_simple(combined_network))

output <- file_path_sans_ext(input_clique_file)
output <- paste0(output, "_enrich.tsv")

network_df <- read_tsv(network_edge_list_file, col_names = FALSE)
ref_genes <- unique(unlist(network_df))
corum_db <- read_gmt(corum_gmt_file)
bioplex_db <- read_gmt(bioplex_gmt_file)

clique_db <- readLines(input_clique_file)
clique_db <- strsplit(clique_db, "\t")

all_results <- NULL

combined_network_vertices <- V(combined_network)$name

for (i in seq_len(length(clique_db))) {
    print(paste0("Processing clique ", i, " of ", length(clique_db)))
    clique_id <- i

    # for each clique, find out what percentage of the edges
    # are covered by any of the interaction databases: BioPlex, BioGrid, HI-Union
    clique_genes <- clique_db[[i]]
    clique_genes_ensembl_id <- unname(symbol_to_id[clique_genes])
    clique_size <- length(clique_genes)
    common_id <- intersect(clique_genes_ensembl_id, combined_network_vertices)
    sub_g <- induced_subgraph(combined_network, common_id)
    sub_g_edges <- gsize(sub_g)
    coverage <- sub_g_edges / (clique_size * (clique_size - 1) / 2.0)
    cur_res <- tibble(clique_id = clique_id,
        clique_list = paste(clique_db[[i]], collapse = ";"),
        clique_size = clique_size,
        clique_coverage_perc = coverage
    )

    # CORUM
    enrich_result <- NULL
    tryCatch(
    enrich_result <- WebGestaltR(enrichMethod = "ORA",
      organism = "hsapiens",
      enrichDatabaseFile = corum_gmt_file,
      enrichDatabaseType = "genesymbol",
      interestGene = clique_db[[i]],
      interestGeneType = "genesymbol",
      referenceGene = ref_genes,
      referenceGeneType = "genesymbol",
      minNum = 3,
      sigMethod = sigmethod,
      topThr = cutoff,
      isOutput = FALSE
    ), error = function(e) {
       print(e)
    })

    if (!is.null(enrich_result)) {
        enrich_result <- as_tibble(enrich_result)
        # should only return 1 row
        stopifnot(dim(enrich_result)[1] == 1)
        cur_gs <- enrich_result[1, "geneSet", drop=TRUE]
        cur_ids <- corum_db[[cur_gs]]
        complex_id <- paste(cur_ids, collapse = ";")
        complex_len <- length(cur_ids)

        colnames(enrich_result) <- paste("corum", colnames(enrich_result), sep = "_")
        cur_res <- cbind(
            cur_res,
            enrich_result,
            corum_complex_list = complex_id,
            corum_complex_size = complex_len
        )
    } else{
        # create empty results
        enrich_result <- tibble(
            geneSet=NA, link=NA, size=NA, overlap=NA, expect=NA,
            enrichmentRatio=NA, pValue=1, FDR=1, overlapId=NA, userId=NA
        )
        colnames(enrich_result) <- paste("corum", colnames(enrich_result), sep = "_")
        cur_res <- cbind(
            cur_res,
            enrich_result,
            corum_complex_list = NA,
            corum_complex_size = NA
        )
    }

    # BioPlex
    enrich_result <- NULL
     tryCatch(
         enrich_result <- WebGestaltR(
             enrichMethod = "ORA",
             organism = "hsapiens",
             enrichDatabaseFile = bioplex_gmt_file,
             enrichDatabaseType = "genesymbol",
             interestGene = clique_db[[i]],
             interestGeneType = "genesymbol",
             referenceGene = ref_genes,
             referenceGeneType = "genesymbol",
             minNum = 3,
             sigMethod = sigmethod,
             topThr = cutoff,
             isOutput = FALSE
         ),
         error = function(e) {
             print(e)
         }
     )

    if (!is.null(enrich_result)) {
        enrich_result <- as_tibble(enrich_result)
        # should only return 1 row
        stopifnot(dim(enrich_result)[1] == 1)
        cur_gs <- enrich_result[1, "geneSet", drop=TRUE]
        cur_ids <- bioplex_db[[cur_gs]]
        complex_id <- paste(cur_ids, collapse = ";")
        complex_len <- length(cur_ids)

        colnames(enrich_result) <- paste("bioplex", colnames(enrich_result), sep = "_")
        cur_res <- cbind(
            cur_res,
            enrich_result,
            bioplex_complex_list = complex_id,
            bioplex_complex_size = complex_len
        )
    } else {
        # create empty results
        enrich_result <- tibble(
            geneSet = NA, link = NA, size = NA, overlap = NA, expect = NA,
            enrichmentRatio = NA, pValue = 1, FDR = 1, overlapId = NA, userId = NA
        )
        colnames(enrich_result) <- paste("bioplex", colnames(enrich_result), sep = "_")
        cur_res <- cbind(
            cur_res,
            enrich_result,
            bioplex_complex_list = NA,
            bioplex_complex_size = NA
        )
    }

   # GO BP
   enrich_result <- NULL
   tryCatch(
    enrich_result <- WebGestaltR(enrichMethod = "ORA",
      organism = "hsapiens",
      enrichDatabase = "geneontology_Biological_Process",
      enrichDatabaseType = "genesymbol",
      interestGene = clique_db[[i]],
      interestGeneType = "genesymbol",
      referenceGene = ref_genes,
      referenceGeneType = "genesymbol",
      minNum = 3,
      sigMethod = sigmethod,
      topThr = cutoff,
      isOutput = FALSE
    ), error = function(e) {
      print(e)
    }
  )
   if (!is.null(enrich_result)) {
        enrich_result <- as_tibble(enrich_result)
        # should only return 1 row
        stopifnot(dim(enrich_result)[1] == 1)
        colnames(enrich_result) <- paste("gobp", colnames(enrich_result), sep = "_")
        cur_res <- cbind(cur_res, enrich_result)
   } else{
        # create empty results
        enrich_result <- tibble(
            geneSet = NA, description=NA, link = NA, size = NA, overlap = NA, expect = NA,
            enrichmentRatio = NA, pValue = 1, FDR = 1, overlapId = NA, userId = NA
        )
        colnames(enrich_result) <- paste("gobp", colnames(enrich_result), sep = "_")
        cur_res <- cbind(
            cur_res,
            enrich_result
        )
   }

  # GO MF
  enrich_result <- NULL
  tryCatch(
      enrich_result <- WebGestaltR(
          enrichMethod = "ORA",
          organism = "hsapiens",
          enrichDatabase = "geneontology_Molecular_Function",
          enrichDatabaseType = "genesymbol",
          interestGene = clique_db[[i]],
          interestGeneType = "genesymbol",
          referenceGene = ref_genes,
          referenceGeneType = "genesymbol",
          minNum = 3,
          sigMethod = sigmethod,
          topThr = cutoff,
          isOutput = FALSE
      ),
      error = function(e) {
          print(e)
      }
  )

   if (!is.null(enrich_result)) {
        enrich_result <- as_tibble(enrich_result)
        # should only return 1 row
        stopifnot(dim(enrich_result)[1] == 1)
        colnames(enrich_result) <- paste("gomf", colnames(enrich_result), sep = "_")
        cur_res <- cbind(cur_res, enrich_result)
   } else {
        # create empty results
        enrich_result <- tibble(
            geneSet = NA,  description=NA, link = NA, size = NA, overlap = NA, expect = NA,
            enrichmentRatio = NA, pValue = 1, FDR = 1, overlapId = NA, userId = NA
        )
        colnames(enrich_result) <- paste("gomf", colnames(enrich_result), sep = "_")
        cur_res <- cbind(cur_res, enrich_result)
   }

    # GO CC
  enrich_result <- NULL
  tryCatch(
      enrich_result <- WebGestaltR(
          enrichMethod = "ORA",
          organism = "hsapiens",
          enrichDatabase = "geneontology_Cellular_Component",
          enrichDatabaseType = "genesymbol",
          interestGene = clique_db[[i]],
          interestGeneType = "genesymbol",
          referenceGene = ref_genes,
          referenceGeneType = "genesymbol",
          minNum = 3,
          sigMethod = sigmethod,
          topThr = cutoff,
          isOutput = FALSE
      ),
      error = function(e) {
          print(e)
      }
  )

   if (!is.null(enrich_result)) {
        enrich_result <- as_tibble(enrich_result)
        # should only return 1 row
        stopifnot(dim(enrich_result)[1] == 1)
        colnames(enrich_result) <- paste("gocc", colnames(enrich_result), sep = "_")
        cur_res <- cbind(cur_res, enrich_result)
   } else {
        # create empty results
        enrich_result <- tibble(
            geneSet = NA,  description=NA, link = NA, size = NA, overlap = NA, expect = NA,
            enrichmentRatio = NA, pValue = 1, FDR = 1, overlapId = NA, userId = NA
        )
        colnames(enrich_result) <- paste("gocc", colnames(enrich_result), sep = "_")
        cur_res <- cbind(cur_res, enrich_result)
   }

   all_results <- rbind(all_results, cur_res)
   write_tsv(all_results, output)
}

# categorize each clique based on pvalues
all_results <- all_results %>%
    mutate(clique_class = if_else(corum_FDR < 0.05, 1,
                    if_else(bioplex_FDR < 0.05, 2,
                    if_else((gobp_FDR < 0.05 | gomf_FDR < 0.05) & gocc_FDR < 0.05, 3, 4))))

write_tsv(all_results, output)
