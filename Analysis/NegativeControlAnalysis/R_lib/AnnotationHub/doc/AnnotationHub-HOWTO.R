## ----style, echo = FALSE, results = 'asis', warning=FALSE-----------------------------------------
options(width=100)
suppressPackageStartupMessages({
    ## load here to avoid noise in the body of the vignette
    library(AnnotationHub)
    library(txdbmaker)
    library(Rsamtools)
    library(VariantAnnotation)
})
BiocStyle::markdown()

## ----less-model-org-------------------------------------------------------------------------------
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "OrgDb")
orgdb <- query(ah, c("OrgDb", "maintainer@bioconductor.org"))[[1]]

## ----less-model-org-select------------------------------------------------------------------------
library(AnnotationDbi)
AnnotationDbi::keytypes(orgdb)
AnnotationDbi::columns(orgdb)
egid <- head(keys(orgdb, "ENTREZID"))
AnnotationDbi::select(orgdb, egid, c("SYMBOL", "GENENAME"), "ENTREZID")

## ----eval=FALSE-----------------------------------------------------------------------------------
# url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E001-H3K4me1.broadPeak.gz"
# filename <-  basename(url)
# download.file(url, destfile=filename)
# if (file.exists(filename))
#    data <- import(filename, format="bed")

## ----results='hide'-------------------------------------------------------------------------------
library(AnnotationHub)
ah = AnnotationHub()
epiFiles <- query(ah, "EpigenomeRoadMap")

## -------------------------------------------------------------------------------------------------
epiFiles

## -------------------------------------------------------------------------------------------------
unique(epiFiles$species)
unique(epiFiles$genome)

## -------------------------------------------------------------------------------------------------
table(epiFiles$sourcetype)

## -------------------------------------------------------------------------------------------------
sort(table(epiFiles$description), decreasing=TRUE)

## -------------------------------------------------------------------------------------------------
metadata.tab <- query(ah , c("EpigenomeRoadMap", "Metadata"))
metadata.tab

## ----echo=FALSE, results='hide'-------------------------------------------------------------------
metadata.tab <- ah[["AH41830"]]

## -------------------------------------------------------------------------------------------------
metadata.tab <- ah[["AH41830"]]

## -------------------------------------------------------------------------------------------------
metadata.tab[1:6, 1:5]

## -------------------------------------------------------------------------------------------------
bpChipEpi <- query(ah , c("EpigenomeRoadMap", "broadPeak", "chip", "consolidated"))

## -------------------------------------------------------------------------------------------------
allBigWigFiles <- query(ah, c("EpigenomeRoadMap", "BigWig"))

## -------------------------------------------------------------------------------------------------
seg <- query(ah, c("EpigenomeRoadMap", "segmentations"))

## -------------------------------------------------------------------------------------------------
E126 <- query(ah , c("EpigenomeRoadMap", "E126", "H3K4ME2"))
E126

## ----echo=FALSE, results='hide'-------------------------------------------------------------------
peaks <- E126[['AH29817']]

## -------------------------------------------------------------------------------------------------
peaks <- E126[['AH29817']]
seqinfo(peaks)

## -------------------------------------------------------------------------------------------------
metadata(peaks)
ah[metadata(peaks)$AnnotationHubName]$sourceurl

## ----takifugu-gene-models-------------------------------------------------------------------------
query(ah, c("Takifugu", "release-94"))

## ----takifugu-data--------------------------------------------------------------------------------
gtf <- ah[["AH64858"]]
dna <- ah[["AH66116"]]

head(gtf, 3)
dna
head(seqlevels(dna))

## ----takifugu-seqlengths--------------------------------------------------------------------------
keep <- names(tail(sort(seqlengths(dna)), 25))
gtf_subset <- gtf[seqnames(gtf) %in% keep]

## ----takifugu-txdb--------------------------------------------------------------------------------
library(txdbmaker)         # for makeTxDbFromGRanges
txdb <- makeTxDbFromGRanges(gtf_subset)

## ----takifugu-exons-------------------------------------------------------------------------------
library(Rsamtools)               # for getSeq,FaFile-method
exons <- exons(txdb)
length(exons)
getSeq(dna, exons)

## -------------------------------------------------------------------------------------------------
chainfiles <- query(ah , c("hg38", "hg19", "chainfile"))
chainfiles

## ----echo=FALSE, results='hide'-------------------------------------------------------------------
chain <- chainfiles[['AH14150']]

## -------------------------------------------------------------------------------------------------
chain <- chainfiles[['AH14150']]
chain

## -------------------------------------------------------------------------------------------------
library(rtracklayer)
gr38 <- liftOver(peaks, chain)

## -------------------------------------------------------------------------------------------------
genome(gr38) <- "hg38"
gr38

## ----echo=FALSE, results='hide', message=FALSE----------------------------------------------------
query(ah, c("GRCh38", "dbSNP", "VCF" ))
vcf <- ah[['AH57960']]

## ----message=FALSE--------------------------------------------------------------------------------
variants <- readVcf(vcf, genome="hg19")
variants

## -------------------------------------------------------------------------------------------------
rowRanges(variants)

## -------------------------------------------------------------------------------------------------
seqlevelsStyle(variants) <-seqlevelsStyle(peaks)

## -------------------------------------------------------------------------------------------------
overlap <- findOverlaps(variants, peaks)
overlap

## -------------------------------------------------------------------------------------------------
idx <- subjectHits(overlap) == 3852
overlap[idx]

## -------------------------------------------------------------------------------------------------
peaks[3852]
rowRanges(variants)[queryHits(overlap[idx])]

## -------------------------------------------------------------------------------------------------
sessionInfo()

