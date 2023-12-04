#' @importFrom digest digest
md5sum <- function(seq){
  md5_name <- vapply(X = as.character(seq),
                     FUN = digest,
                     FUN.VALUE = character(1L),
                     algo="md5")
}

#' @importFrom S4Vectors mcols
checkRMSK <- function(rmsk){
  stopifnot(is(rmsk, "GRanges"))
  if(!all(.globals$repColNames %in% colnames(mcols(rmsk)))){
    stop('rmsk must be a GRanges object',
         'with meta data of repName, repClass, and repFamily')
  }
}

#' @importFrom methods is
checkInputs <- function(genome, txdb, rmsk, peaks){
  checkRMSK(rmsk)
  stopifnot(is(genome, 'BSgenome'))
  stopifnot(is(txdb, 'TxDb'))
  if(!missing(peaks)){
    stopifnot(is(peaks, 'GRanges'))
  }
}

#' @importFrom data.table .N .SD as.data.table
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom IRanges findOverlaps
prepareAnno <- function(seq_gl, rmsk, seq,
                        seq_type=.globals$feature_types){
  checkRMSK(rmsk)
  stopifnot(is(seq_gl, 'GRangesList'))
  seq_type <- match.arg(seq_type)
  ol <- findOverlaps(query = rmsk, subject = seq_gl)
  rm_anno <- rmsk[queryHits(ol)]
  rm_anno$feature <- seq_type
  rm_anno$seqname <- names(seq)[subjectHits(ol)]
  rm_anno$idx <- as.character(rm_anno)
  if(length(names(seq_gl))==length(seq_gl)){
    rm_anno$tx_name <- names(seq_gl)[subjectHits(ol)]
  }else{
    rm_anno$tx_name <- NA
  }
  rm_anno <- mcols(rm_anno)[
    c('seqname', 'tx_name', 'feature', .globals$repColNames, 'idx')]
  rm_anno <- as.data.table(rm_anno)
  rm_anno <- rm_anno[, c(.N, lapply(.SD, paste, collapse=',')),
                     by=c('seqname', 'feature',
                          .globals$repColNames),
                     .SDcols = c('idx', 'tx_name')]
  rm_anno <- DataFrame(rm_anno)
  rm_anno$seq <- seq[rm_anno$seqname]
  rm_anno
}

#' @importFrom methods as
rmSeqByName <- function(seq, gr, seqnames){
  stopifnot(is(seq, 'DNAStringSet'))
  stopifnot(is(gr, 'GRanges'))
  k <- names(seq) %in% seqnames
  if(any(k)){
    gr <- gr[-k]
    seq <- seq[-k]
  }
  if(!is.null(gr$tx_name) && length(gr$tx_name)==length(gr)){
    n <- vapply(gr$tx_name, `[`, i=1, FUN.VALUE = character(1L))
    gr <- as(gr, 'GRangesList')
    names(gr) <- n
  }else{
    gr <- as(gr, 'GRangesList')
  }
  list(gr=gr, seq=seq)
}

# private function, input will be checked in upstream function
prepareSeqWithAnno <- function(rmsk, genome, gr,
                               seq_type, seqnamesToBeRemoved,
                               subsetGRanges,
                               ...){
  if(!missing(subsetGRanges)){
    gr <- subsetByOverlaps(gr, subsetGRanges)
  }
  seq <- getSeq(x = genome, names = gr)
  names(seq) <- md5sum(seq)
  fil <- rmSeqByName(seq, gr, seqnamesToBeRemoved)
  anno <- prepareAnno(fil$gr, rmsk, fil$seq, seq_type = seq_type)
  anno
}

#' @importFrom GenomicFeatures exonsBy extractTranscriptSeqs
#' @importFrom IRanges subsetByOverlaps
getExonsSeq <- function(genome, txdb, rmsk, subsetGRanges, ...){
  checkInputs(genome, txdb, rmsk)
  exons <- exonsBy(x = txdb, by = 'tx', use.names=TRUE)
  if(!missing(subsetGRanges)){
    exons <- subsetByOverlaps(exons, subsetGRanges)
  }
  exonSeq <- extractTranscriptSeqs(x=genome, transcripts=exons)
  names(exonSeq) <- md5sum(exonSeq)
  anno <- prepareAnno(exons, rmsk, exonSeq, seq_type = 'exon')
  anno
}

#' @importFrom GenomicFeatures intronicParts
#' @importFrom BSgenome getSeq
getIntronsSeq <- function(genome, txdb, rmsk, seqnamesToBeRemoved, ...){
  checkInputs(genome, txdb, rmsk)
  if(missing(seqnamesToBeRemoved)){
    stop('seqnamesToBeRemoved is required')
  }
  introns <- intronicParts(txdb = txdb) # this set the priority of exon first
  prepareSeqWithAnno(rmsk, genome, introns, seq_type = 'intron',
                     seqnamesToBeRemoved, ...)
}

reMinWidth <- function(gr, minWidth, txdb){
  gr_width <- width(gr)
  k <- gr_width < minWidth
  center_ <- start(gr[k]) + k[k]/2
  st <- floor(center_ - minWidth/2)
  start(gr[k]) <- ifelse(st<1, 1, st)
  se <- ceiling(center_ + minWidth/2)
  l <- seqlengths(x = txdb)[as.character(seqnames(gr[k]))]
  end(gr[k]) <- ifelse(se>l, l, se)
  gr
}

#' @importFrom GenomicFeatures genes
#' @importMethodsFrom IRanges subsetByOverlaps reduce
#' @importFrom BiocGenerics start end width `start<-` `end<-`
#' @importFrom GenomeInfoDb seqnames seqlengths
getIntergenicSeq <- function(genome, txdb, rmsk,
                             seqnamesToBeRemoved,
                             minWidth=31,
                             ...){
  checkInputs(genome, txdb, rmsk)
  if(missing(seqnamesToBeRemoved)){
    stop('seqnamesToBeRemoved is required')
  }
  genes <- genes(x = txdb, single.strand.genes.only = FALSE)
  intergenic <- subsetByOverlaps(x = rmsk, ranges = genes, invert = TRUE)
  intergenic <- reduce(x = intergenic, ignore.strand = TRUE)
  intergenic <- reMinWidth(intergenic, minWidth=minWidth, txdb = txdb)
  prepareSeqWithAnno(rmsk, genome, intergenic, seq_type = 'intergenic',
                     seqnamesToBeRemoved, ...)
}

#' @importFrom IRanges nearest
#' @importFrom GenomicFeatures transcripts
getPeakSeq <- function(genome, txdb, rmsk, peaks,
                       minWidth=31,
                       ...){
  checkInputs(genome, txdb, rmsk, peaks)
  peaks <- reduce(peaks)
  txs <- transcripts(txdb)
  if(length(txs$tx_name)==length(txs)){
    nid <- nearest(peaks, txs)
    peaks$tx_name <- txs$tx_name[nid]
  }
  peaks <- reMinWidth(peaks, minWidth=minWidth, txdb = txdb)
  prepareSeqWithAnno(rmsk, genome, peaks, seq_type = 'peak',
                     seqnamesToBeRemoved=c(), ...)
}

isPrepareSeqOut <- function(seq){
  msg <- 'Input is not a output from prepareSeq function.'
  if(!is(seq, 'DataFrame')){
    stop(msg)
  }
  if(!'seq' %in% colnames(seq)){
    stop(msg)
  }
}

