#' rollup single salmon quants
#' @noRd
#' @param quant_sf The salmon quant.sf file path.
#' @param seq_anno Output of \link{prepareSeq} function.
#' @return A list of counts table with columns TPM and NumReads
#' @importFrom readr read_tsv
#' @importFrom dplyr `%>%` group_by_at summarize_at left_join
#' @importFrom tibble column_to_rownames
rollup_one_RE <- function(quant_sf, seq_anno){
  stopifnot(is.character(quant_sf))
  stopifnot(file.exists(quant_sf))
  qsf <- read_tsv(quant_sf, col_types="cdddd") # c: "character", d: "double"
  stopifnot(
    "Unexpected colnames of salmon quants file."=
      all(.globals$quantSfColName %in% colnames(qsf)))
  qsf <- qsf[, .globals$quantSfColName]
  seqCnts <- as.data.frame(seq_anno) %>%
    left_join(qsf, by=c("seqname"="Name"))
  repCnts <- lapply(.globals$repColNamesRollup, FUN = function(.col){
    seqCnts %>%
      group_by_at(.vars=.col) %>%
      summarize_at(.globals$quantSfColName[-1], sum, na.rm=TRUE) %>%
      column_to_rownames(var=.col)
  })
  repCnts
}
#' Rollup multiple salmon quants
#' @description
#' Rollup the TPM and NumReads for repeat elements from salmon quants.
#' 
#' @param quant_sfs The salmon quant.sf files path.
#' @param qnames Human-readable names for each quant.sf file.
#' @param seq_anno Output of \link{prepareSeq} function.
#' @param select_feature Type of feature. Available choices are:
#' 'all', 'exon', 'intron', 'gene', or 'intergenic'.
#' @param cpus Number of threads.
#' @param retry Retry number for the rollup.
#' @importFrom parallel mclapply
#' @export
#' @examples
#' quant_sfs <- dir(system.file('extdata', 'RNA-seq', package='ratte'),
#'         'quant.sf',
#'  full.names=TRUE)
#' qnames <- sub("^.*\\/(.*?).quant.sf", "\\1", quant_sfs)
#' rmsk <- readRDS(system.file('extdata', 'danRer11_rmsk_chr25_sub.RDS',
#'                              package='ratte'))
#' 
#' library(BSgenome.Drerio.UCSC.danRer11) # for sequences
#' library(TxDb.Drerio.UCSC.danRer11.refGene) # for genomic annotation
#' seq_anno <- prepareSeq(Drerio,
#'            TxDb.Drerio.UCSC.danRer11.refGene,
#'            rmsk,
#'            subsetGRanges=GRanges('chr25:11000000-15000000'))
#' re <- rollup_RE(quant_sfs, qnames, seq_anno, 'all')
#' re_gene_level <- rollup_RE(quant_sfs, qnames, seq_anno, 'gene')
rollup_RE <- function(quant_sfs, qnames, seq_anno,
                      select_feature,
                      cpus=2,
                      retry=3){
  stopifnot("quant_sfs is empty."=length(quant_sfs)>=1)
  stopifnot(length(quant_sfs)==length(qnames))
  stopifnot(is.numeric(retry))
  stopifnot(retry>0)
  names(quant_sfs) <- qnames
  select_feature <- 
    match.arg(select_feature,
              choices = .globals$feature_types_choices)
  seq_anno <- switch (select_feature,
    'all' = seq_anno,
    'exon' = seq_anno[seq_anno$feature=='exon', , drop=FALSE],
    'intron' = seq_anno[seq_anno$feature=='intron', , drop=FALSE],
    'intergenic' = seq_anno[seq_anno$feature=='intergenic', , drop=FALSE],
    'gene' = seq_anno[seq_anno$feature %in% c('exon', 'intron'), , drop=FALSE],
    'peak' = seq_anno[seq_anno$feature=='peak', , drop=FALSE]
  )
  see <- TRUE
  while(retry>0 && any(see)){
    retry <- retry - 1
    se <- mclapply(quant_sfs, FUN = function(quant_sf){
      rollup_one_RE(quant_sf = quant_sf, seq_anno = seq_anno)
    }, mc.cores=cpus)
    see <- mapply(inherits, se, 'try-error')
  }
  if(any(see)){
    message(se[see])
    stop('Error in read the quant files. You may want to retry it again.')
  }
  rn <- rownames(se[[1]][[1]])
  rn_identical <- vapply(se, FUN = function(.ele){
    all(rownames(.ele[[1]])==rn)
  }, FUN.VALUE = logical(1L))
  stopifnot("Rownames of quant file are not identical."=all(rn_identical))
  coln <- .globals$quantSfColName[-1]
  se <- lapply(coln, FUN = function(.col){
    lapply(.globals$repColNamesRollup, FUN = function(.repN){
      x <- do.call(cbind, lapply(se, FUN = function(.ele){
        .ele[[.repN]][, .col, drop=FALSE]
      }))
      colnames(x) <- qnames
      x
    })
  })
  se <- unlist(se, recursive = FALSE)
}

#' Create a list of edgeR DGEList object from salmon quants
#' @description
#' Create a list of edgeR DGEList object for repeat elements from salmon quants.
#' 
#' @param quant_sfs The salmon quant.sf files path.
#' @param qnames Human-readable names for each quant.sf file.
#' @param seq_anno Output of \link{prepareSeq} function.
#' @param select_feature Type of feature. Available choices are:
#' 'all', 'exon', 'intron', 'gene', or 'intergenic'.
#' @param cpus Number of threads.
#' @param norm_method Normalization method for edgeR.
#' @param norm_level Normalization feature level for edgeR.
#' @param ... Parameters will be passed to \link[edgeR:DGEList]{DGEList}.
#' Please note, the default norm.factors and lib.size will be calculated
#' at gene level.
#' @return A list of DGEList objects.
#' @importFrom edgeR DGEList calcNormFactors
#' @export
#' @examples
#' quant_sfs <- dir(system.file('extdata', 'RNA-seq', package='ratte'),
#'         'quant.sf',
#'  full.names=TRUE)
#' qnames <- sub("^.*\\/(.*?).quant.sf", "\\1", quant_sfs)
#' rmsk <- readRDS(system.file('extdata', 'danRer11_rmsk_chr25_sub.RDS',
#'                              package='ratte'))
#' 
#' library(BSgenome.Drerio.UCSC.danRer11) # for sequences
#' library(TxDb.Drerio.UCSC.danRer11.refGene) # for genomic annotation
#' seq_anno <- prepareSeq(Drerio,
#'            TxDb.Drerio.UCSC.danRer11.refGene,
#'            rmsk,
#'            subsetGRanges=GRanges('chr25:11000000-15000000'))
#' ys <- createDGELists(quant_sfs, qnames, seq_anno, 'all',
#'                      group=sub('.rep.', '', qnames))
#' ys[[1]]
#' names(ys)
#' library(edgeR)
#' res <- lapply(ys, function(y){
#'   keep <- filterByExpr(y)
#'   y <- y[keep, keep.lib.sizes=TRUE]
#'   design <- model.matrix(~y$samples$group)
#'   y <- estimateDisp(y, design)
#'   fit <- glmFit(y, design)
#'   lrt <- glmLRT(fit, coef=2)
#'   topTags(lrt, n=nrow(lrt))
#' })
#' tt <- as.data.frame(res[[1]])
#' plot(tt$logFC, -10*log10(tt$PValue),
#'    xlab='logFC', ylab='-10logP')
createDGELists <- function(quant_sfs, qnames, seq_anno,
                           select_feature, cpus=2,
                           norm_method = "RLE",
                           norm_level = "gene",
                           ...){
  norm_level <- 
    match.arg(norm_level,
              choices = .globals$feature_types_choices)
  re <- rollup_RE(quant_sfs, qnames, seq_anno, select_feature, cpus)
  if(!norm_level %in% seq_anno$feature &&
     !(norm_level=='gene' && any(c('exon', 'intron') %in% seq_anno$feature))){
    warning(norm_level, ' is not available in the annotation. Use "all".')
    norm_level <- 'all'
  }
  if(select_feature!=norm_level){
    re_gene <- rollup_RE(quant_sfs, qnames, seq_anno, norm_level, cpus)
  }else{
    re_gene <- re
  }
  y <- lapply(.globals$repColNamesRollup, function(.col){
    .n <- paste(.globals$quantSfColName[3], .col, sep=".")
    dots <- list(...)
    dots$counts <- re[[.n]]
    if(is.null(dots$libSize)){
      dots$libSize <- colSums(re_gene[[.n]], na.rm=TRUE)
    }
    if(is.null(dots$normFactors)){
      dots$normFactors <- calcNormFactors(object = re_gene[[.n]],
                                          method = norm_method)
    }
    if(.col=='seqname'){
      dots$gene <- seq_anno[match(rownames(dots$counts),
                                  seq_anno$seqname), , drop=FALSE]
      rownames(dots$gene) <- dots$gene[, 1]
      dots$gene <- dots$gene[, -1, drop=FALSE]
    }
    dgel <- do.call(DGEList, dots)
    dgel$tpm <- re[[paste(.globals$quantSfColName[2],
                          .col, sep=".")]][rownames(dgel), ]
    dgel
  })
  names(y)[names(y)=='seqname'] <- 'rawCounts'
  y
}
