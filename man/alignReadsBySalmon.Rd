% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareAnnotations.R
\name{alignReadsBySalmon}
\alias{alignReadsBySalmon}
\title{Quatifying the samples by salmon}
\usage{
alignReadsBySalmon(
  R1,
  R2,
  output,
  salmonPath = "salmon",
  cpus = 2,
  index = "TEindex",
  dryrun = TRUE,
  ...
)
}
\arguments{
\item{R1, R2}{The R1 and/or R2 reads file path.}

\item{output}{The output folder for salmon.}

\item{salmonPath}{The path to salmon executable file.}

\item{cpus}{Threads number.}

\item{index}{Output index folder.}

\item{dryrun}{If true, print the command only.}

\item{...}{Parameters passed to system2.}
}
\value{
If dryrun is false, the output of system2.
}
\description{
Run salmon quant command for the given fastqs.
}
\examples{
alignReadsBySalmon('sample1.R1.fq.gz',
                   output='sample1.salmon.REdiscoverTE',
                   dryrun=TRUE)

}
