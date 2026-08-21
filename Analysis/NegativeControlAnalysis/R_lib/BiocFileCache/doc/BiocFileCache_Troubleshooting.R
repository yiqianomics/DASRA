## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ----eval = FALSE-------------------------------------------------------------
# if (!"BiocManager" %in% rownames(installed.packages()))
#      install.packages("BiocManager")
# BiocManager::install("BiocFileCache", dependencies=TRUE)

## ----results='hide', warning=FALSE, message=FALSE-----------------------------
library(BiocFileCache)

## ----create-------------------------------------------------------------------
path <- tempfile()
bfc <- BiocFileCache(path, ask = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#        # make sure you have permissions on the cache/files
#        # use at own risk
# 
# 	moveFiles<-function(package){
# 	    olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
# 	    newdir <- tools::R_user_dir(package, which="cache")
# 	    dir.create(path=newdir, recursive=TRUE)
# 	    files <- list.files(olddir, full.names =TRUE)
# 	    moveres <- vapply(files,
# 		FUN=function(fl){
# 		  filename = basename(fl)
# 		  newname = file.path(newdir, filename)
# 		  file.rename(fl, newname)
# 		},
# 		FUN.VALUE = logical(1))
# 	    if(all(moveres)) unlink(olddir, recursive=TRUE)
# 	}
# 
# 
# 	package="BiocFileCache"
# 	moveFiles(package)
# 

## ----eval=FALSE---------------------------------------------------------------
# library(BiocFileCache)
# 
# 
# package = "BiocFileCache"
# 
# BFC_CACHE = rappdirs::user_cache_dir(appname=package)
# Sys.setenv(BFC_CACHE = BFC_CACHE)
# bfc = BiocFileCache(BFC_CACHE)
# ## CAUTION: This removes the cache and all downloaded resources
# removebfc(bfc, ask=FALSE)
# 
# ## create new empty cache in new default location
# bfc = BiocFileCache(ask=FALSE)
# 

## ----sessioninfo--------------------------------------------------------------
sessionInfo()


