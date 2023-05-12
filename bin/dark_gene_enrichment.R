#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
library(argparse)
library(tools)
library(httr)
library(stats)
library(igraph)
library(jsonlite)
library(parallel)
suppressPackageStartupMessages(library(WebGestaltR))


netwalker <- function(seed, adjMatrix, r = 0.5) {
  de <- apply(adjMatrix, 2, sum)
  w <- t(t(adjMatrix)/de)

  p0 <- array(0, dim = c(nrow(w), 1))
  rownames(p0) <- rownames(w)
  p0[seed, 1] <- 1/length(seed)

  pt <- p0
  pt1 <- (1-r)*(w%*%pt)+r*p0

  while (sum(abs(pt1-pt)) > 1e-6) {
    pt <- pt1
    pt1 <- (1-r)*(w%*%pt)+r*p0
  }

  return(pt1)
}


get_top_neighbors <- function(seeds, id, net_adj, net, top_n, prefix) {
  results <- list()
  for (i in seq_along(seeds)) {
    cat(paste0(id, ": [", i, "/", length(seeds), "] ",
      ": ", seeds[i], "\n"))
    seed <- c(seeds[i])
    pt1 <- netwalker(seed, net_adj, r = 0.5)
    net_node <- V(net)$name
    gS <- data.frame(name = net_node, score = pt1, stringsAsFactors = F)
    gS_filtered <- gS[!(gS$name %in% list(seed)), ]
    gS_sorted <- gS_filtered[order(gS_filtered$score, decreasing = TRUE), ]
    candidate <- gS_sorted[1:top_n, ]
    selected <- candidate$name
    results[[seed]] <- selected
    rm(gS, gS_sorted, gS_filtered, candidate)
  }

  df <- do.call(cbind.data.frame, results)
  write.table(df, file = paste0(prefix, "_", id, ".tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  return(df)
}


dark_enrich <- function(id, top_neighor_prefix, enrich_results_prefix) {
    all_results <- NULL
    top_neighor_file <- paste0(top_neighor_prefix, "_", id, ".tsv")
    # read from file
    df <- read.table(top_neighor_file, header = FALSE, sep = "\t",
        stringsAsFactors = FALSE)
    n_seeds <- ncol(df)
    for (i in seq_len(n_seeds)) {
        print(paste0(i, " of ", n_seeds))
        enrich_result <- NULL
        interest_gene <- df %>% pull(i)
        dark_gene <- interest_gene[1]
        tryCatch(
        enrich_result <- WebGestaltR(enrichMethod = "ORA",
            organism = "hsapiens",
            enrichDatabase = "geneontology_Biological_Process",
            enrichDatabaseType = "genesymbol",
            interestGene = interest_gene,
            interestGeneType = "genesymbol",
            referenceGene = net_node,
            referenceGeneType = "genesymbol",
            minNum = 10,
            sigMethod = "top",
            tppThr = 10,
            isOutput = FALSE
        ), error = function(e) {
        print(e)
        }
        )
        if (!is.null(enrich_result)) {
        cur_res <- cbind(
                dark_gene = dark_gene,
                dark_gene_top_neighbor = paste(interest_gene, collapse = ";"),
                top_neighbor_size = length(interest_gene),
                ontology = "BP",
                enrich_result
        )
        all_results <- rbind(all_results, cur_res)
        }

        enrich_result <- NULL
        tryCatch(
        enrich_result <- WebGestaltR(enrichMethod = "ORA",
            organism = "hsapiens",
            enrichDatabase = "geneontology_Cellular_Component",
            enrichDatabaseType = "genesymbol",
            interestGene = interest_gene,
            interestGeneType = "genesymbol",
            referenceGene = net_node,
            referenceGeneType = "genesymbol",
            minNum = 10,
            sigMethod = "top",
            topThr = 10,
            isOutput = FALSE
        ), error = function(e) {
        print(e)
        }
        )
        if (!is.null(enrich_result)) {
        cur_res <- cbind(
                dark_gene = dark_gene,
                dark_gene_top_neighbor = paste(interest_gene, collapse = ";"),
                top_neighbor_size = length(interest_gene),
                ontology = "CC",
                enrich_result
        )
        all_results <- rbind(all_results, cur_res)
        }

        enrich_result <- NULL
        tryCatch(
        enrich_result <- WebGestaltR(enrichMethod = "ORA",
            organism = "hsapiens",
            enrichDatabase = "geneontology_Molecular_Function",
            enrichDatabaseType = "genesymbol",
            interestGene = interest_gene,
            interestGeneType = "genesymbol",
            referenceGene = net_node,
            referenceGeneType = "genesymbol",
            minNum = 10,
            sigMethod = "top",
            topThr = 10,
            isOutput = FALSE
        ), error = function(e) {
        print(e)
        }
        )
        if (!is.null(enrich_result)) {
        cur_res <- cbind(
                dark_gene = dark_gene,
                dark_gene_top_neighbor = paste(interest_gene, collapse = ";"),
                top_neighbor_size = length(interest_gene),
                ontology = "MF",
                enrich_result
        )
        all_results <- rbind(all_results, cur_res)
        }
  }
  cur_file <- paste0(enrich_results_prefix, "_", id, ".tsv")
  write.table(all_results, file = cur_file,
    sep = "\t", quote = FALSE, row.names = FALSE)
  return(cur_file)
}


## main
parser <- ArgumentParser()
parser$add_argument("--cut-off-count", type = "integer",
    default = 0,
    help = "the maximum number of publication to be considered as dark genes"
)
parser$add_argument("--n-cpus", type = "integer",
    default = 1,
    help = "number of cpus available"
)
parser$add_argument("--top-n-neighbor", type = "integer",
  default = 50,
  help = paste0("for the darker gene, how many top scored NTA to",
           "consider for the enrichment analysis")
)
parser$add_argument("--network-el", type = "character",
  required = TRUE,
  help = paste0("network edgelist file for NTA analysis")
)
parser$add_argument("--gene-pubmed-count", type = "character",
  help = paste0("path the gene pubmed count file")
)


args <- parser$parse_args()
cut_off_count <- args$cut_off_count
n_cpus <- args$n_cpus
top_n_neighbor <- args$top_n_neighbor
network_el <- args$network_el
start <- args$start
end <- args$end
gene_pubmed <- args$gene_pubmed_count

df <- read_tsv(network_el, col_names = FALSE, col_types = "cc")
all_genes_in_network <- unique(unlist(df))

g2p_df <- read_tsv(gene_pubmed, col_types = "ci")
colnames(g2p_df) <- c("gene_symbol", "count")
g2p_df <- g2p_df %>% drop_na()

cnf_df <- tibble(gene_symbol = all_genes_in_network)
cnt_df <- cnf_df %>% left_join(g2p_df, by = "gene_symbol")
cnt_df <- cnt_df %>% replace_na(list(count = 0))

# * what is qualified as dark gene?  with publication count <= cut_off_count
dark_genes <- cnt_df %>%
            filter(count <= cut_off_count) %>%
            pull("gene_symbol")
dark_genes <- sort(dark_genes)
print(length(dark_genes))

# find the neighbors of the dark genes and write to files
net <- as.matrix(read_tsv(network_el, col_names = FALSE, col_types = "cc"))
net_graph <- graph.edgelist(net, directed = FALSE)
# all dark genes must already be in the network
# sanity check
net_node <- V(net_graph)$name
inter <- intersect(net_node, dark_genes)
if (length(dark_genes) != length(inter)) {
stop(paste0("The number of dark genes is not equal to the number",
        " of nodes in the network!"))
}

cl <- makeCluster(n_cpus, outfile = "")
#   Calculate the number of genes per chunk
genes_per_chunk <- ceiling(length(dark_genes) / n_cpus)
gene_chunks <- split(dark_genes, rep(1:n_cpus, each = genes_per_chunk,
    length.out = length(dark_genes)))
n_chunks <- length(gene_chunks)
top_neighor_prefix <- "top_neighbor"

ids <- seq_len(n_chunks)
# make it 2 digits
ids <- sprintf("%02d", ids)
net_adjmatrix <- get.adjacency(net_graph, sparse = FALSE)
clusterExport(cl, c("netwalker", "get.adjacency", "V", "net_adjmatrix"))
arg_list <- list(data = gene_chunks, net_adj = net_adjmatrix, net = net_graph,
     top_n = top_n_neighbor, prefix = top_neighor_prefix, id = ids)
results <- mclapply(seq_along(arg_list$data),
        function(i)  get_top_neighbors(
        arg_list$data[[i]],
        arg_list$id[i],
        arg_list$net_adj,
        arg_list$net,
        arg_list$top_n,
        arg_list$prefix
        ), mc.cores = length(cl))

stopCluster(cl)
# concatenate the files into one column wise
# list all files  in the current directory that matches the pattern
# top_neighbor_*.tsv
top_neighbor_files <- list.files(pattern = paste0(top_neighor_prefix, "_.*tsv"))
top_neighbor_files <- sort(top_neighbor_files)
# use system command to concatenate the files
cmd <- paste0("paste ", paste(top_neighbor_files, collapse = " "), " > ",
    top_neighor_prefix, ".tsv")
system(cmd)


# do the enrichment analysis
cl <- makeCluster(n_cpus, outfile = "")
clusterExport(cl, c("V", "%>%", "pull", "WebGestaltR"))

enrich_results_prefix <- "enrich_results"
arg_list <- list(
     id = ids, top_neighor_prefix = top_neighor_prefix,
        enrich_results_prefix = enrich_results_prefix
)
results <- mclapply(seq_along(ids),
        function(i)  dark_enrich(
        arg_list$id[i],
        arg_list$top_neighor_prefix,
        arg_list$enrich_results_prefix
        ), mc.cores = length(cl))

stopCluster(cl)

# combine the files into one with name enrich_results_*.tsv
# list all files  in the current directory that matches the pattern
# enrich_results_*.tsv
enrich_files <- list.files(pattern = paste0(enrich_results_prefix, "_.*tsv"))
# sort the files
enrich_files <- sort(enrich_files)
# read all files and combine them into one file
all_enrich_results <- lapply(enrich_files, read_tsv)
# write to file
write.table(do.call(rbind, all_enrich_results), file = "enrich_results.tsv",
    sep = "\t", quote = FALSE, row.names = FALSE)
