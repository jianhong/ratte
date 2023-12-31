% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareAnnotations.R
\name{prepareSeq}
\alias{prepareSeq}
\title{Prepare the sequences with repeat elements}
\usage{
prepareSeq(genome, txdb, rmsk, peaks, minWidth = 31L, verbose = FALSE, ...)
}
\arguments{
\item{genome}{A BSgenome object}

\item{txdb}{A TxDb object}

\item{rmsk}{The RepeatMasker GRanges object}

\item{peaks}{An GRanges object to represent the peak list.}

\item{minWidth}{Minimal width for the repeat elements in the intergenic
annotation.}

\item{verbose}{Print the message or not.}

\item{...}{Not use.}
}
\value{
An DataFrame with seqname and their repeat element annotations
}
\description{
The function will extract the genomic annotations from TxDb object,
extract the sequences from BSgenome object, and
generate the annotation for repeat elements.
}
\examples{
# library(AnnotationHub)
# ah <- AnnotationHub()
# (rmsk <- query(ah, c('RepeatMasker', 'Zebrafish')))
# rmsk <- rmsk[['AH98980']]
rmsk <- readRDS(system.file('extdata', 'danRer11_rmsk_chr25_sub.RDS',
                             package='ratte'))

library(BSgenome.Drerio.UCSC.danRer11) # for sequences
library(TxDb.Drerio.UCSC.danRer11.refGene) # for genomic annotation
prepareSeq(Drerio,
           TxDb.Drerio.UCSC.danRer11.refGene,
           rmsk)
}
