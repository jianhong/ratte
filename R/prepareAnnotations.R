#' Prepare the sequences with repeat elements
#' @description
#' The function will extract the genomic annotations from TxDb object,
#' extract the sequences from BSgenome object, and
#' generate the annotation for repeat elements.
#' @param genome A BSgenome object
#' @param txdb A TxDb object
#' @param rmsk The RepeatMasker GRanges object
#' @param peaks An GRanges object to represent the peak list.
#' @param minWidth Minimal width for the repeat elements in the intergenic
#'  annotation.
#' @param verbose Print the message or not.
#' @param ... Not use.
#' @return An DataFrame with seqname and their repeat element annotations
#' @export
#' @examples
#' # library(AnnotationHub)
#' # ah <- AnnotationHub()
#' # (rmsk <- query(ah, c('RepeatMasker', 'Zebrafish')))
#' # rmsk <- rmsk[['AH98980']]
#' rmsk <- readRDS(system.file('extdata', 'danRer11_rmsk_chr25_sub.RDS',
#'                              package='ratte'))
#' 
#' library(BSgenome.Drerio.UCSC.danRer11) # for sequences
#' library(TxDb.Drerio.UCSC.danRer11.refGene) # for genomic annotation
#' prepareSeq(Drerio,
#'            TxDb.Drerio.UCSC.danRer11.refGene,
#'            rmsk)
prepareSeq <- function(genome, txdb, rmsk, peaks,
                       minWidth=31L, verbose = FALSE, ...){
  checkInputs(genome, txdb, rmsk, peaks)
  stopifnot(is.numeric(minWidth))
  stopifnot(is.logical(verbose))
  if(!missing(peaks)){
    if(verbose) message('Prepare the peaks ...')
    seq <- getPeakSeq(genome, txdb, rmsk, peaks, minWidth = minWidth, ...)
  }else{
    if(verbose) message('Prepare the exons...')
    exonSeq <- getExonsSeq(genome, txdb, rmsk, ...)
    if(verbose) message('Prepare the introns...')
    intronSeq <- getIntronsSeq(genome, txdb, rmsk, exonSeq$seqname, ...)
    if(verbose) message('Prepare the intergenics...')
    intergenicSeq <- getIntergenicSeq(genome, txdb, rmsk,
                                      c(exonSeq$seqname, intronSeq$seqname),
                                      minWidth = minWidth,
                                      ...)
    seq <- rbind(exonSeq, intronSeq, intergenicSeq)
  }
  seq
}

#' Export the sequences as fasta file.
#' @description
#' Export the sequences as a gzip fasta file.
#' @param seq The output of \link{prepareSeq} function.
#' @param filepath The output fasta file path.
#' @return Logical value to indicate the success of output.
#' @importFrom Biostrings writeXStringSet
#' @export
#' @examples
#' # library(AnnotationHub)
#' # ah <- AnnotationHub()
#' # (rmsk <- query(ah, c('RepeatMasker', 'Zebrafish')))
#' # rmsk <- rmsk[['AH98980']]
#' rmsk <- readRDS(system.file('extdata', 'danRer11_rmsk_chr25_sub.RDS',
#'                              package='ratte'))
#' library(BSgenome.Drerio.UCSC.danRer11) # for sequences
#' library(TxDb.Drerio.UCSC.danRer11.refGene) # for genomic annotation
#' seq <- prepareSeq(Drerio,
#'            TxDb.Drerio.UCSC.danRer11.refGene,
#'            rmsk,
#'            subsetGRanges=GRanges('chr25:11000000-15000000'))
#' outfile <- tempfile(fileext='.fa.gz')
#' saveFasta(seq, outfile)
saveFasta <- function(seq, filepath='rmsk.fa.gz'){
  isPrepareSeqOut(seq)
  if(!grepl('(fa|fasta).gz', filepath)){
    stop('filepath must end with fa.gz or fasta.gz')
  }
  writeXStringSet(x = unique(seq$seq), filepath = filepath, format='fasta',
                  compress = 'gzip', compression_level = 9)
}

#' Create the index files for salmon
#' @description
#' Run salmon index command for the given fasta.
#' @param salmonPath The path to salmon executable file.
#' @param cpus Threads number.
#' @param FaFile The fasta file.
#' @param index Output index folder.
#' @param dryrun If true, print the command only.
#' @param ... Parameters passed to system2.
#' @return If dryrun is false, the output of system2.
#' @export
#' @examples
#' prepareSalmonIndex('salmon', dryrun=TRUE)
#' 
prepareSalmonIndex <- function(salmonPath='salmon',
                               cpus=2,
                               FaFile='rmsk.fa.gz',
                               index='TEindex',
                               dryrun=TRUE,
                               ...){
  stopifnot(is.numeric(cpus))
  args <- c('index', '--threads', cpus, '-t', FaFile, '-i', index)
  if(dryrun){
    cmd <- paste(args, collapse = " ")
    message(salmonPath, " ", cmd)
  }else{
    stopifnot(file.exists(FaFile))
    system2(command=salmonPath,
            args = args,
            ...)
  }
}

#' Quatifying the samples by salmon
#' @description
#' Run salmon quant command for the given fastqs.
#' @param R1,R2 The R1 and/or R2 reads file path.
#' @param output The output folder for salmon.
#' @param salmonPath The path to salmon executable file.
#' @param cpus Threads number.
#' @param index Output index folder.
#' @param dryrun If true, print the command only.
#' @param ... Parameters passed to system2.
#' @return If dryrun is false, the output of system2.
#' @export
#' @examples
#' alignReadsBySalmon('sample1.R1.fq.gz',
#'                    output='sample1.salmon.REdiscoverTE',
#'                    dryrun=TRUE)
#' 
alignReadsBySalmon <- function(R1, R2, output,
                               salmonPath='salmon',
                               cpus=2,
                               index='TEindex',
                               dryrun=TRUE,
                               ...){
  stopifnot(is.numeric(cpus))
  if(missing(R1) || missing(output)){
    stop('R1 and output is requred.')
  }
  args <- c('quant', '--seqBias', '--gcBias', 
            '--index', index,
            '--libType', 'A',
            '--validateMappings',
            '--threads', cpus,
            '-o', output)
  if(missing(R2)){
    args <- c(args, '--unmatedReads', R1)
  }else{
    args <- c(args, '-1', R1, '-2', R2)
  }
  if(dryrun){
    cmd <- paste(args, collapse = " ")
    message(salmonPath, " ", cmd)
  }else{
    stopifnot(file.exists(R1))
    stopifnot(file.exists(R2))
    system2(command=salmonPath,
            args = args,
            ...)
  }
}

#' Prepare repeat element annotations
#' @description
#' Create unique annotation table for repeat element to their class and 
#' families.
#' @param rmsk The RepeatMasker GRanges object
#' @return An DataFrame object
#' @export
#' @examples
#' rmsk <- readRDS(system.file('extdata', 'danRer11_rmsk_chr25_sub.RDS',
#'                              package='ratte'))
#' prepareREannotation(rmsk)
prepareREannotation <- function(rmsk){
  checkRMSK(rmsk)
  unique(mcols(rmsk)[, .globals$repColNames])
}
