# ratte
The use of RnA-seq for AnnoTation of Transposable Elements (ratte) is a package
to enhance the power of [REdiscoverTE](https://github.com/ucsffrancislab/REdiscoverTE)
by break down the barrier of species. 
To broaden the scope of annotation to encompass more species, the ratte package 
leverages the resources available within the Bioconductor community to prepare
the necessary sequence files and annotations.

# Installation

```{r}
library(devtools)
devtools::install("jianhong/ratte")
```

# Online documentation

The full documentation is available [here](https:://jianhong.github.io/ratte).

```{r}
## load libraries
library(ratter)
library(AnnotationHub)
library(BSgenome.Drerio.UCSC.danRer11)
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(edgeR)

## get RepeatMasker elements
ah <- AnnotationHub()
(rmsk <- query(ah, c('RepeatMasker', 'Zebrafish')))
rmsk <- rmsk[['AH98980']]

## prepare the annotation and fasta file
seq_anno <- prepareSeq(Drerio, TxDb.Drerio.UCSC.danRer11.refGene, rmsk)
saveFasta(seq_anno, 'danRer11.rmsk.fasta.gz')

## create the salmon index
prepareSalmonIndex(salmonPath='salmon',
                   FaFile='danRer11.rmsk.fasta.gz',
                   index='TEindex')
## align the reads
fastqs <- list(R1=c('WT1_R1.fq.gz', 'WT2_R1.fq.gz',
                    'KD1_R1.fq.gz', 'KD2_R1.fq.gz'),
               R2=c('WT1_R2.fq.gz', 'WT2_R2.fq.gz',
                    'KD1_R2.fq.gz', 'KD2_R2.fq.gz'),
               N=c('WT1', 'WT2', 'KD1', 'KD2'))
output_path <- 'salmon_output'
dir.create(output_path)
for(i in seq_along(fastqs[[1]])){
alignReadsBySalmon(R1=fastqs[['R1']][i],
                   R2=fastqs[['R2']][i],
                   output=file.path(output_path, fastqs[['N']][i]),
                   index='TEindex')
}

## Aggregate alignments to repeat element
quant_sfs <- dir(output_path, 'quant.sf',
                 full.names=TRUE,
                 recursive=TRUE)
qnames <- c('WT1', 'WT2', 'KD1', 'KD2')
group <- c('WT', 'WT', 'KD', 'KD')
ys <- createDGELists(quant_sfs, qnames, seq_anno,
                     select_feature = 'all',
                     group=group,
                     remove.zeros = TRUE)

## Differential analysis
y <- ys[['repName']]
keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes=TRUE]
design <- model.matrix(~y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
res <- as.data.frame(topTags(lrt, n=nrow(lrt)))
```


