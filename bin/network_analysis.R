#!/usr/bin/env Rscript

library(NetSAM)
library(argparse)
library(tools)
library(GO.db)
library(DBI)
library(tidyverse)
library(igraph)

# modified from NetSAM::GOAssociation, output data, not html
GO_Association <- function(NetSAMOutput, outputFile,
    organism="hsapiens",
    outputType="significant", fdrmethod="BH", fdrth=0.05,
    topNum=5){

    if (is.null(NetSAMOutput$rulfile) || is.null(NetSAMOutput$hmifile)) {
        stop("NetSAMOutput should be the output from NetSAM function, which contain rulfile, hmifile and network!\n")
    }

    if(length(which(outputType %in% c("significant","top")))==0){
        stop("The input 'outputType' is invalid! Please select a method from 'significant' and 'top'!\n")
    }

    organisms <- c("hsapiens","mmusculus","rnorvegicus","drerio","celegans","scerevisiae","cfamiliaris","dmelanogaster","athaliana")
    names(organisms) <- c("Hs","Mm","Rn","Dr","Ce","Sc","Cf","Dm","At")

    if(length(which(organisms==organism))==0){
        stop("Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana!")
    }

    or <- names(organisms)[which(organisms==organism)]

    if(length(which(fdrmethod %in% c("holm","hochberg","hommel","bonferroni","BH","BY","none")))==0){
        stop("The input 'fdrmethod' is invalid! Please select an option from 'holm','hochberg', 'hommel', 'bonferroni', 'BH', 'BY' and 'none'!\n")
    }

    if(organism=="scerevisiae"){
        .sql_bp <-  "select distinct t1.gene_name,t2.go_id from sgd as t1 inner join go_bp_all as t2 on t1._id=t2._id"
        .sql_cc <-  "select distinct t1.gene_name,t2.go_id from sgd as t1 inner join go_cc_all as t2 on t1._id=t2._id"
        .sql_mf <-  "select distinct t1.gene_name,t2.go_id from sgd as t1 inner join go_mf_all as t2 on t1._id=t2._id"
    }else{
        .sql_bp <-  "select distinct t1.symbol,t2.go_id from gene_info as t1 inner join go_bp_all as t2 on t1._id=t2._id"
        .sql_cc <-  "select distinct t1.symbol,t2.go_id from gene_info as t1 inner join go_cc_all as t2 on t1._id=t2._id"
        .sql_mf <-  "select distinct t1.symbol,t2.go_id from gene_info as t1 inner join go_mf_all as t2 on t1._id=t2._id"
    }
    if(organism=="scerevisiae"){
        conn <- org.Sc.sgd.db::org.Sc.sgd_dbconn()
    }else{
        if(organism=="athaliana"){
            conn <- org.At.tair.db::org.At.tair_dbconn()
        }else{
            conn <- eval(parse(text=paste("org.",or,".eg.db::","org.",or,".eg_dbconn()", sep = "")))
        }
    }
    generalAnn_BP <- dbGetQuery(conn, .sql_bp)
    generalAnn_BP <- generalAnn_BP[!is.na(generalAnn_BP[,1]),]
    bp_num <- tapply(generalAnn_BP[,1],generalAnn_BP[,2],length)
    bp_num <- bp_num[bp_num>=10 & bp_num<=2000]
    generalAnn_BP <- generalAnn_BP[generalAnn_BP[,2] %in% names(bp_num),]

    generalAnn_CC <- dbGetQuery(conn, .sql_cc)
    generalAnn_CC <- generalAnn_CC[!is.na(generalAnn_CC[,1]),]
    cc_num <- tapply(generalAnn_CC[,1],generalAnn_CC[,2],length)
    cc_num <- cc_num[cc_num>=10 & cc_num<=2000]
    generalAnn_CC <- generalAnn_CC[generalAnn_CC[,2] %in% names(cc_num),]

    generalAnn_MF <- dbGetQuery(conn, .sql_mf)
    generalAnn_MF <- generalAnn_MF[!is.na(generalAnn_MF[,1]),]
    mf_num <- tapply(generalAnn_MF[,1],generalAnn_MF[,2],length)
    mf_num <- mf_num[mf_num>=10 & mf_num<=2000]
    generalAnn_MF <- generalAnn_MF[generalAnn_MF[,2] %in% names(mf_num),]


    .sql <- "select distinct go_id goid,term name from go_term";
    conn <- get("GO_dbconn")()
    allTermName <- dbGetQuery(conn,.sql)

    rul <- NetSAMOutput$rulfile
    hmi <- NetSAMOutput$hmifile

    refGenes <- rul[,4]

    annRef_BP <- generalAnn_BP[generalAnn_BP[,1] %in% refGenes,]
    annRef_CC <- generalAnn_CC[generalAnn_CC[,1] %in% refGenes,]
    annRef_MF <- generalAnn_MF[generalAnn_MF[,1] %in% refGenes,]

    moduleEnrich <- data.frame(module_name = "", module_size = 0,
      module_overlap = 0,
      ontology = "", gene_set = "", gene_set_desc = "",
      gene_set_size = 0,
      ref_overlap = 0,
      pvalue = 0, FDR = 0, stringsAsFactors = F)
    mi <- 1
    for(i in c(2:nrow(hmi))){
        mN <- hmi[i, 4]
        s <- hmi[i, 5]
        e <- hmi[i, 6]
        sz <- e - s + 1
        g <- rul[s:e, 4]

        annInterest_BP <- generalAnn_BP[generalAnn_BP[, 1] %in% g, ]
        annInterest_CC <- generalAnn_CC[generalAnn_CC[, 1] %in% g, ]
        annInterest_MF <- generalAnn_MF[generalAnn_MF[, 1] %in% g, ]

        t <- c()
        if(nrow(annInterest_BP)>0){
            termInfo <- .enrichmentFunction(annRef_BP, annInterest_BP,allTermName,fdrmethod)
            if(outputType=="significant"){
                termInfo <- termInfo[termInfo[,6]<fdrth,]
                if(nrow(termInfo)>0){
                    termInfo$ontology <- "BP"
                    termInfo$go_size <- bp_num[termInfo$goid]
                    t <- rbind(t,termInfo)
                }
            }else{
                termInfo <- termInfo[order(termInfo[,5]),]
                termInfo <- termInfo[1:topNum,]
                termInfo$ontology <- "BP"
                termInfo$go_size <- bp_num[termInfo$goid]
                t <- rbind(t,termInfo)
            }
        }
        if(nrow(annInterest_CC)>0){
            termInfo <- .enrichmentFunction(annRef_CC, annInterest_CC,allTermName,fdrmethod)
            if(outputType=="significant"){
                termInfo <- termInfo[termInfo[,6]<fdrth,]
                if(nrow(termInfo)>0){
                    termInfo$ontology <- "CC"
                    termInfo$go_size <- cc_num[termInfo$goid]
                    t <- rbind(t,termInfo)
                }
            }else{
                termInfo <- termInfo[order(termInfo[,5]),]
                termInfo <- termInfo[1:topNum,]
                termInfo$ontology <- "CC"
                termInfo$go_size <- cc_num[termInfo$goid]
                t <- rbind(t,termInfo)
            }
        }
        if(nrow(annInterest_MF)>0){
            termInfo <- .enrichmentFunction(annRef_MF, annInterest_MF,
              allTermName,fdrmethod)
            if(outputType=="significant"){
                termInfo <- termInfo[termInfo[,6]<fdrth,]
                if(nrow(termInfo)>0){
                    termInfo$ontology <- "MF"
                    termInfo$go_size <- mf_num[termInfo$goid]
                    t <- rbind(t,termInfo)
                }
            }else{
                termInfo <- termInfo[order(termInfo[,5]),]
                termInfo <- termInfo[1:topNum,]
                termInfo$go_size <- mf_num[termInfo$goid]
                termInfo$ontology <- "MF"
                t <- rbind(t,termInfo)
            }
        }

        if(length(t)>0){
            moduleEnrich[mi:(mi+nrow(t)-1),1] <- mN
            moduleEnrich[mi:(mi+nrow(t)-1),2] <- sz
            moduleEnrich[mi:(mi+nrow(t)-1),3] <- t[,4]
            moduleEnrich[mi:(mi+nrow(t)-1),4] <- t[,7]
            moduleEnrich[mi:(mi+nrow(t)-1),5] <- t[,1]
            moduleEnrich[mi:(mi+nrow(t)-1),6] <- t[,2]
            moduleEnrich[mi:(mi+nrow(t)-1),7] <- t[,8]
            moduleEnrich[mi:(mi+nrow(t)-1),8] <- t[,3]
            moduleEnrich[mi:(mi+nrow(t)-1),9] <- format(t[,5],
              scientific=TRUE,digits=3)
            moduleEnrich[mi:(mi+nrow(t)-1),10] <- format(t[,6],
              scientific=TRUE,digits=3)
            mi <- mi + nrow(t)
        }
    }

    if (moduleEnrich[1, 1] != "") {
       write.table(moduleEnrich, file = outputFile, sep = "\t",
          row.names = FALSE, col.names = TRUE, quote = FALSE)
        return(moduleEnrich)
    } else{
      cat(paste0("There is no associated GO term ",
          "for each module based on the input parameters!\n"))
      return(NULL)
    }

}


.enrichmentFunction <- function(annRef, annInterest, allTermName, fdrmethod)
{
    allRefnum <- length(unique(annRef[,1]))
    allInterestnum <- length(unique(annInterest[,1]))

    allAnnterm <- unique(annRef[,2])
    allAnntermL <- length(allAnnterm)


    refTermCount <- tapply(annRef[,1],annRef[,2],length)

    refTermName <- allTermName[allTermName[,1] %in% unique(annRef[,2]),]
    rownames(refTermName) <- refTermName[,1]
    refTermName <- refTermName[names(refTermCount),]

    refTermCount <- data.frame(goid=names(refTermCount),name=refTermName[,2],refnum=refTermCount,stringsAsFactors=F)
    refTermCount <- refTermCount[order(refTermCount[,1]),]
    interestTermCount <- tapply(annInterest[,1],annInterest[,2],length)
    interestTermCount <- data.frame(goid=names(interestTermCount),interestnum=interestTermCount,stringsAsFactors=F)
    interestTermCount <- interestTermCount[order(interestTermCount[,1]),]

    ref_interest_TermCount <- refTermCount

    ref_interest_TermCount$interestnum = array(0,dim=c(length(ref_interest_TermCount$goid),1))
    ref_interest_TermCount[ref_interest_TermCount$goid %in% interestTermCount[,1],4]=interestTermCount$interestnum


    n <- nrow(ref_interest_TermCount)
    pv <- array(0,dim=c(n,1))
    for (i in c(1:n)){
        p <- 1-phyper(ref_interest_TermCount[i,4]-1,allInterestnum,allRefnum-allInterestnum,ref_interest_TermCount[i,3],lower.tail = TRUE,log.p= FALSE)
        pv[i,1] <- p
    }
    ref_interest_TermCount$pvalue <- pv
    adp <- p.adjust(pv,method=fdrmethod)
    ref_interest_TermCount$FDR <- adp
    return(ref_interest_TermCount)
}




## main function
parser <- ArgumentParser()
parser$add_argument("--input", type = "character", required = TRUE,
    help = "path to the network edge list file",
    metavar = "path_to_input_network")
parser$add_argument("--netsam-output-file", type = "character",
    help = "(optional) previoulsly generated netsam output file")
parser$add_argument("--min-module-size", type = "integer", default = 5,
    help = "minimum module size (default: 5)")
parser$add_argument("--enrich-top-n",
    type = "integer", default = 10,
    help = "number of top enriched terms to output (default: 10)")
parser$add_argument("--nthreads",
    type = "integer", default = 4,
    help = "number of threads (default: 4)")

args <- parser$parse_args()
input_network <- args$input
min_module_size <- args$min_module_size
top_n <- args$enrich_top_n
nthreads <- as.numeric(args$nthreads)

network <- read.graph(input_network, format = "ncol")
min_module_ratio <- min_module_size / vcount(network)

# check if netsam output is already generated
if (is.null(args$netsam_output_file)) {
  dir_name <- dirname(input_network)
#   current_date_time <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
  out_dir <- file.path(dir_name, "out")
  dir.create(out_dir, showWarnings = FALSE)
  output_filename <- paste0(file_path_sans_ext(basename(input_network)),
      "_out_", min_module_size)
  output_file <- file.path(out_dir, output_filename)
  netsam_out_file <- paste0(output_file, ".RData")
} else{
  netsam_out_file <- args$netsam_output_file
  output_file <- file.path(dirname(nestam_out_file),
      file_path_sans_ext(basename(nestam_out_file)))
}


if (is.null(args$netsam_output_file)) {
  netsam_out <- NetSAM(inputNetwork = input_network, outputFileName = output_file,
    outputFormat = "nsm", edgeType = "unweighted", map_to_genesymbol = FALSE,
    organism = "hsapiens", idType = "auto", minModule = min_module_ratio,
    stepIte = TRUE, maxStep = 10, moduleSigMethod = "permutation",
    modularityThr = 0.2,
    ZRanNum = 10, PerRanNum = 1000, ranSig = 0.05, edgeThr = (-1),
    nodeThr = (-1), nThreads = nthreads)
  save(netsam_out, file = netsam_out_file)
  NetAnalyzer(input_network, output_file, "unweighted")
} else {
    cat("loading", netsam_out_file, "...")
    load(netsam_out_file)
}

# performa go association analysis for all cases
# top_n <- 1
# go_enrich_output_file <- paste0(output_file, "_go_enrich_top_", top_n, ".tsv")
# GO_Association(NetSAMOutput = netsam_out, outputType = "top",
#   outputFile = go_enrich_output_file,
#   organism = "hsapiens", topNum = top_n)
go_enrich_output_file <- paste0(output_file, "_go_enrich_top_", top_n, ".tsv")
GO_Association(NetSAMOutput = netsam_out, outputType = "top",
  outputFile = go_enrich_output_file,
  organism = "hsapiens", topNum = top_n)
