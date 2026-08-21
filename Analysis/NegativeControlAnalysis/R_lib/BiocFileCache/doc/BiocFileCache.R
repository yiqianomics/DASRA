## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ----eval = FALSE-------------------------------------------------------------
# if (!"BiocManager" %in% rownames(installed.packages()))
#      install.packages("BiocManager")
# BiocManager::install("BiocFileCache", dependencies=TRUE)

## ----library, results='hide', warning=FALSE, message=FALSE--------------------
library(BiocFileCache)

## ----create-------------------------------------------------------------------
path <- tempfile()
bfc <- BiocFileCache(path, ask = FALSE)

## ----cacheloc-----------------------------------------------------------------
bfccache(bfc)
length(bfc)

## ----bfcshow------------------------------------------------------------------
bfc

## ----bfcinfo------------------------------------------------------------------
bfcinfo(bfc)

## ----bfcnew-------------------------------------------------------------------
savepath <- bfcnew(bfc, "NewResource", ext=".RData")
savepath

## now we can use that path in any save function
m = matrix(1:12, nrow=3)
save(m, file=savepath)

## and that file will be tracked in the cache
bfcinfo(bfc)

## ----bfcadd-------------------------------------------------------------------
fl1 <- tempfile(); file.create(fl1)
add2 <- bfcadd(bfc, "Test_addCopy", fl1)                 # copy
# returns filepath being tracked in cache
add2
# the name is the unique rid in the cache
rid2 <- names(add2)

fl2 <- tempfile(); file.create(fl2)
add3 <- bfcadd(bfc, "Test2_addMove", fl2, action="move") # move
rid3 <- names(add3)

fl3 <- tempfile(); file.create(fl3)
add4 <- bfcadd(bfc, "Test3_addAsis", fl3, rtype="local",
	       action="asis") # reference
rid4 <- names(add4)

file.exists(fl1)    # TRUE - copied from original location
file.exists(fl2)    # FALSE - moved from original location
file.exists(fl3)    # TRUE - left asis, original location tracked

## ----bfcaddremote-------------------------------------------------------------
url <- "http://httpbin.org/get"
add5 <- bfcadd(bfc, "TestWeb", fpath=url)
rid5 <- names(add5)

url2<- "https://bioconductor.org/packages/stats/bioc/BiocFileCache/BiocFileCache_2024_stats.tab"
add6 <- bfcadd(bfc, "TestWeb", fpath=url2)
rid6 <- names(add6)

# add a remote resource but don't initially download
add7 <- bfcadd(bfc, "TestNoDweb", fpath=url2, download=FALSE)
rid7 <- names(add7)
# let's look at our BiocFileCache object now
bfc
bfcinfo(bfc)

## ----bfcquery-----------------------------------------------------------------
bfcquery(bfc, "Web")

bfcquery(bfc, "copy")

q1 <- bfcquery(bfc, "BiocFileCache")
q1
class(q1)

## ----bfccount-----------------------------------------------------------------
bfccount(q1)

## ----bfcsubset----------------------------------------------------------------
bfcsubWeb = bfc[paste0("BFC", 5:6)]
bfcsubWeb
bfcinfo(bfcsubWeb)

## ----bfcbracket---------------------------------------------------------------
bfc[["BFC2"]]
bfcpath(bfc, "BFC2")
bfcpath(bfc, "BFC5")
bfcrpath(bfc, rids="BFC5")
bfcrpath(bfc)
bfcrpath(bfc, c("http://httpbin.org/get","Test3_addAsis"))

## ----bfcneedsupdate-----------------------------------------------------------
bfcneedsupdate(bfc, "BFC5")
bfcneedsupdate(bfc, "BFC6")
bfcneedsupdate(bfc)

## ----bfcrename----------------------------------------------------------------
fileBeingReplaced <- bfc[[rid3]]
fileBeingReplaced

# fl3 was created when we were adding resources
fl3

bfc[[rid3]]<-fl3
bfc[[rid3]]

## ----bfcupdate----------------------------------------------------------------
bfcinfo(bfc, "BFC1")
bfcupdate(bfc, "BFC1", rname="FirstEntry")
bfcinfo(bfc, "BFC1")

## ----bfcupdateremote----------------------------------------------------------
suppressPackageStartupMessages({
    library(dplyr)
})
bfcinfo(bfc, "BFC6") %>% select(rid, rpath, fpath)
bfcupdate(bfc, "BFC6", fpath=url, rname="Duplicate", ask=FALSE)
bfcinfo(bfc, "BFC6") %>% select(rid, rpath, fpath)

## ----bfcdownload--------------------------------------------------------------
rid <- "BFC5"
test <- !identical(bfcneedsupdate(bfc, rid), FALSE) # 'TRUE' or 'NA'
if (test)
    bfcdownload(bfc, rid, ask=FALSE)

## ----bfcmetadata--------------------------------------------------------------
names(bfcinfo(bfc))
meta <- as.data.frame(list(rid=bfcrid(bfc)[1:3], idx=1:3))
bfcmeta(bfc, name="resourceData") <- meta
names(bfcinfo(bfc))

## ----bfcmetalist--------------------------------------------------------------
bfcmetalist(bfc)
bfcmeta(bfc, name="resourceData")

## ----bfcmetaremove------------------------------------------------------------
bfcmetaremove(bfc, name="resourceData")

## ----eval=FALSE---------------------------------------------------------------
# bfcmeta(name="resourceData") <- meta
# Error in bfcmeta(name = "resourceData") <- meta :
#   target of assignment expands to non-language object

## ----eval=FALSE---------------------------------------------------------------
# bfc <- BiocFileCache()
# bfcmeta(bfc, name="resourceData") <- meta

## ----bfcremove----------------------------------------------------------------
# let's remind ourselves of our object
bfc

bfcremove(bfc, "BFC6")
bfcremove(bfc, "BFC1")

# let's look at our BiocFileCache object now
bfc

## ----bfcsync------------------------------------------------------------------
# create a new entry that hasn't been used
path <- bfcnew(bfc, "UseMe")
rmMe <- names(path)
# We also have a file not being tracked because we updated rpath

bfcsync(bfc)

# you can suppress the messages and just have a TRUE/FALSE
bfcsync(bfc, FALSE)

#
# Let's do some cleaning to have a synced object
#
bfcremove(bfc, rmMe)
unlink(fileBeingReplaced)

bfcsync(bfc)

## ----eval=FALSE---------------------------------------------------------------
# # export entire biocfilecache
# exportbfc(bfc)
# 
# # export the first 4 entries of biocfilecache
# # as a compressed tar
# exportbfc(bfc, rids=paste0("BFC", 1:4),
# 	  outputFile="BiocFileCacheExport.tar.gz", compression="gzip")
# 
# # export the subsetted object of web resources as zip
# sub1 <- bfc[bfcrid(bfcquery(bfc, "web", field='rtype'))]
# exportbfc(sub1, outputFile = "BiocFileCacheExportWeb.zip",
# 	  outMethod="zip")

## ----eval=FALSE---------------------------------------------------------------
# 
# bfc <- importbfc("BiocFileCacheExport.tar")
# 
# bfc2 <- importbfc("BiocFileCacheExport.tar.gz", compression="gzip")
# 
# bfc3 <- importbfc("BiocFileCacheExportWeb.zip", archiveMethod="unzip")

## ----mock---------------------------------------------------------------------
tbl <- data.frame(rtype=c("web","web"),
		      rpath=c(NA_character_,NA_character_),
		  fpath=c("http://httpbin.org/get",
			  "https://en.wikipedia.org/wiki/Bioconductor"),
		      keywords = c("httpbin", "wiki"), stringsAsFactors=FALSE)
tbl

## ----eval=FALSE---------------------------------------------------------------
# 
# newbfc <- makeBiocFileCacheFromDataFrame(tbl,
# 					 cache=file.path(tempdir(),"BFC"),
# 					 actionWeb="copy",
# 					 actionLocal="copy",
# 					 metadataName="resourceMetadata")
# 

## ----eval=FALSE---------------------------------------------------------------
# cleanbfc(bfc)

## ----eval=FALSE---------------------------------------------------------------
# removebfc(bfc)

## ----eval=FALSE---------------------------------------------------------------
# proxy <- httr::use_proxy("http://my_user:my_password@myproxy:8080")
# ## or
# proxy <- httr::use_proxy(Sys.getenv('http_proxy'))
# httr::set_config(proxy)

## ----sessioninfo--------------------------------------------------------------
sessionInfo()


