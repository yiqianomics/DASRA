## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ----eval = FALSE-------------------------------------------------------------
# if (!"BiocManager" %in% rownames(installed.packages()))
#      install.packages("BiocManager")
# BiocManager::install("BiocFileCache", dependencies=TRUE)

## ----results='hide', warning=FALSE, message=FALSE-----------------------------
library(BiocFileCache)

## -----------------------------------------------------------------------------
path <- tempfile()
bfc <- BiocFileCache(path, ask = FALSE)

## ----url----------------------------------------------------------------------
## paste to avoid long line in vignette
url <- paste(
    "ftp://ftp.ensembl.org/pub/release-71/gtf",
    "homo_sapiens/Homo_sapiens.GRCh37.71.gtf.gz",
    sep="/")

## ----eval=FALSE---------------------------------------------------------------
# library(BiocFileCache)
# bfc <- BiocFileCache()
# path <- bfcrpath(bfc, url)

## ----eval=FALSE---------------------------------------------------------------
# gtf <- rtracklayer::import.gff(path)

## ----eval=FALSE---------------------------------------------------------------
# gtf <- rtracklayer::import.gff(bfcrpath(BiocFileCache(), url))

## ----eval=FALSE---------------------------------------------------------------
# library(BiocFileCache)
# bfc <- BiocFileCache("~/my-experiment/results")

## ----eval=FALSE---------------------------------------------------------------
# suppressPackageStartupMessages({
#     library(DESeq2)
#     library(airway)
# })
# data(airway)
# dds <- DESeqDataData(airway, design = ~ cell + dex)
# result <- DESeq(dds)

## ----eval=FALSE---------------------------------------------------------------
# saveRDS(result, bfcnew(bfc, "airway / DESeq standard analysis"))

## ----eval=FALSE---------------------------------------------------------------
# result <- readRDS(bfcrpath(bfc, "airway / DESeq standard analysis"))

## ----eval=FALSE---------------------------------------------------------------
# suppressPackageStartupMessages({
#     library(BiocFileCache)
#     library(rtracklayer)
# })
# 
# # load the cache
# path <- file.path(tempdir(), "tempCacheDir")
# bfc <- BiocFileCache(path)
# 
# # the web resource of interest
# url <- "ftp://ftp.ensembl.org/pub/release-71/gtf/homo_sapiens/Homo_sapiens.GRCh37.71.gtf.gz"
# 
# # check if url is being tracked
# res <- bfcquery(bfc, url, exact=TRUE)
# 
# if (bfccount(res) == 0L) {
# 
#     # if it is not in cache, add
#     ans <- bfcadd(bfc, rname="ensembl, homo sapien", fpath=url)
# 
# } else {
# 
#   # if it is in cache, get path to load
#   rid = res$rid
#   ans <- bfcrpath(bfc, rid)
# 
#   # check to see if the resource needs to be updated
#   check <- bfcneedsupdate(bfc, rid)
#   # check can be NA if it cannot be determined, choose how to handle
#   if (is.na(check)) check <- TRUE
#   if (check){
#     ans < - bfcdownload(bfc, rid)
#   }
# }
# 
# # ans is the path of the file to load
# ans
# 
# # we know because we search for the url that the file is a .gtf.gz,
# # if we searched on other terms we can use 'bfcpath' to see the
# # original fpath to know the appropriate load/read/import method
# bfcpath(bfc, names(ans))
# 
# temp = GTFFile(ans)
# info = import(temp)

## ----ensemblremote, eval=TRUE-------------------------------------------------

#
# A simpler test to see if something is in the cache
# and if not start tracking it is using `bfcrpath`
#

suppressPackageStartupMessages({
    library(BiocFileCache)
    library(rtracklayer)
})

# load the cache
path <- file.path(tempdir(), "tempCacheDir")
bfc <- BiocFileCache(path, ask=FALSE)

# the web resources of interest
url <- "ftp://ftp.ensembl.org/pub/release-71/gtf/homo_sapiens/Homo_sapiens.GRCh37.71.gtf.gz"

url2 <- "ftp://ftp.ensembl.org/pub/release-71/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_5.0.71.gtf.gz"

# if not in cache will download and create new entry
pathsToLoad <- bfcrpath(bfc, c(url, url2))

pathsToLoad

# now load files as see fit
info = import(GTFFile(pathsToLoad[1]))
class(info)
summary(info)

## ----eval=FALSE---------------------------------------------------------------
# #
# # One could also imagine the following:
# #
# 
# library(BiocFileCache)
# 
# # load the cache
# bfc <- BiocFileCache()
# 
# #
# # Do some work!
# #
# 
# # add a location in the cache
# filepath <- bfcnew(bfc, "R workspace")
# 
# save(list = ls(), file=filepath)
# 
# # now the R workspace is being tracked in the cache

## ----eval=FALSE---------------------------------------------------------------
# .get_cache <-
#     function()
# {
#     cache <- tools::R_user_dir("MyNewPackage", which="cache")
#     BiocFileCache::BiocFileCache(cache)
# }

## ----eval=FALSE---------------------------------------------------------------
# download_data_file <-
#     function( verbose = FALSE )
# {
#     fileURL <- "http://a_path_to/someremotefile.tsv.gz"
# 
#     bfc <- .get_cache()
#     rid <- bfcquery(bfc, "geneFileV2", "rname")$rid
#     if (!length(rid)) {
# 	 if( verbose )
# 	     message( "Downloading GENE file" )
# 	 rid <- names(bfcadd(bfc, "geneFileV2", fileURL ))
#     }
#     if (!isFALSE(bfcneedsupdate(bfc, rid)))
# 	bfcdownload(bfc, rid)
# 
#     bfcrpath(bfc, rids = rid)
# }

## ----preprocess---------------------------------------------------------------
url <- "http://bioconductor.org/packages/stats/bioc/BiocFileCache/BiocFileCache_stats.tab"

headFile <-                         # how to process file before caching
    function(from, to)
{
    dat <- readLines(from)
    writeLines(head(dat), to)
    TRUE
}

rid <- bfcquery(bfc, url, "fpath")$rid
if (!length(rid))                   # not in cache, add but do not download
    rid <- names(bfcadd(bfc, url, download = FALSE))

update <- bfcneedsupdate(bfc, rid)  # TRUE if newly added or stale
if (!isFALSE(update))               # download & process
    bfcdownload(bfc, rid, ask = FALSE, FUN = headFile)

rpath <- bfcrpath(bfc, rids=rid)    # path to processed result
readLines(rpath)                    # read processed result

## ----sessioninfo--------------------------------------------------------------
sessionInfo()


