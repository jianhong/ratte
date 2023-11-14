.globals <- new.env(parent = emptyenv())
.globals$repColNames <- c('repName', 'repClass', 'repFamily')
names(.globals$repColNames) <- .globals$repColNames
.globals$quantSfColName <- c('Name', 'TPM', 'NumReads')
names(.globals$quantSfColName) <- .globals$quantSfColName
