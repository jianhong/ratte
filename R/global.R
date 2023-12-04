.globals <- new.env(parent = emptyenv())
.globals$repColNames <- c('repName', 'repClass', 'repFamily')
names(.globals$repColNames) <- .globals$repColNames
.globals$repColNamesRollup <- c(.globals$repColNames, 'seqname'='seqname')
.globals$quantSfColName <- c('Name', 'TPM', 'NumReads')
names(.globals$quantSfColName) <- .globals$quantSfColName
.globals$feature_types_choices <-
  c('all', 'exon', 'intron', 'gene', 'intergenic', 'peak')
.globals$feature_types <- 
  c('exon', 'intron', 'intergenic', 'peak')