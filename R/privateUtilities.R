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
                        seq_type=c('exon', 'intron', 'intergenic')){
  checkRMSK(rmsk)
  stopifnot(is(seq_gl, 'GRangesList'))
  seq_type <- match.arg(seq_type)
  ol <- findOverlaps(query = rmsk, subject = seq_gl)
  rm_anno <- rmsk[queryHits(ol)]
  rm_anno$feature <- seq_type
  rm_anno$seqname <- names(seq)[subjectHits(ol)]
  rm_anno$idx <- as.character(rm_anno)
  rm_anno <- mcols(rm_anno)[
    c('seqname', 'feature', .globals$repColNames, 'idx')]
  rm_anno <- as.data.table(rm_anno)
  rm_anno <- rm_anno[, c(.N, lapply(.SD, paste, collapse=',')),
                     by=c('seqname', 'feature',
                          .globals$repColNames),
                     .SDcols = c('idx')]
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
  list(gr=as(gr, 'GRangesList'), seq=seq)
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
getExonsSeq <- function(genome, txdb, rmsk, peaks, subsetGRanges, ...){
  checkInputs(genome, txdb, rmsk, peaks)
  if(!missing(peaks)){
    exons <- exons(x = txdb)
    exons <- subsetByOverlaps(x = exons, ranges = peaks, minoverlap = 1L,
                              ignore.strand = TRUE)
  }else{
    exons <- exonsBy(x = txdb, by = 'tx', use.names=TRUE)
  }
  if(!missing(subsetGRanges)){
    exons <- subsetByOverlaps(exons, subsetGRanges)
  }
  if(!missing(peaks)){
    exonSeq <- getSeq(x=genome, names=exons)
    exons <- as(exons, 'GRangesList')
  }else{
    exonSeq <- extractTranscriptSeqs(x=genome, transcripts=exons)
  }
  names(exonSeq) <- md5sum(exonSeq)
  anno <- prepareAnno(exons, rmsk, exonSeq, seq_type = 'exon')
  anno
}

#' @importFrom GenomicFeatures intronicParts
#' @importFrom BSgenome getSeq
getIntronsSeq <- function(genome, txdb, rmsk, peaks, seqnamesToBeRemoved, ...){
  checkInputs(genome, txdb, rmsk, peaks)
  if(missing(seqnamesToBeRemoved)){
    stop('seqnamesToBeRemoved is required')
  }
  introns <- intronicParts(txdb = txdb) # this set the priority of exon first
  if(!missing(peaks)){
    introns <- subsetByOverlaps(x = introns, ranges = peaks, minoverlap = 1L,
                                ignore.strand = TRUE)
  }
  prepareSeqWithAnno(rmsk, genome, introns, seq_type = 'intron',
                     seqnamesToBeRemoved, ...)
}

#' @importFrom GenomicFeatures genes
#' @importMethodsFrom IRanges subsetByOverlaps reduce
#' @importFrom BiocGenerics start end width `start<-` `end<-`
#' @importFrom GenomeInfoDb seqnames seqlengths
getIntergenicSeq <- function(genome, txdb, rmsk, peaks,
                             seqnamesToBeRemoved,
                             minWidth=31,
                             ...){
  checkInputs(genome, txdb, rmsk, peaks)
  if(missing(seqnamesToBeRemoved)){
    stop('seqnamesToBeRemoved is required')
  }
  genes <- genes(x = txdb, single.strand.genes.only = FALSE)
  intergenic <- subsetByOverlaps(x = rmsk, ranges = genes, invert = TRUE)
  intergenic <- reduce(x = intergenic, ignore.strand = TRUE)
  if(!missing(peaks)){
    intergenic <- subsetByOverlaps(x = intergenic, ranges = peaks,
                                   minoverlap = 1L,
                                   ignore.strand = TRUE)
  }
  intergenic_width <- width(intergenic)
  k <- intergenic_width < minWidth
  center_ <- start(intergenic[k]) + k[k]/2
  st <- floor(center_ - minWidth/2)
  start(intergenic[k]) <- ifelse(st<1, 1, st)
  se <- ceiling(center_ + minWidth/2)
  l <- seqlengths(x = txdb)[as.character(seqnames(intergenic[k]))]
  end(intergenic[k]) <- ifelse(se>l, l, se)
  prepareSeqWithAnno(rmsk, genome, intergenic, seq_type = 'intergenic',
                     seqnamesToBeRemoved, ...)
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

