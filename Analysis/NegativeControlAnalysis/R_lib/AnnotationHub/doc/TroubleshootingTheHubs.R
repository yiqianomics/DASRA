## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(eval = FALSE)

## ----eval=FALSE-----------------------------------------------------------------------------------
# 
#        # make sure you have permissions on the cache/files
#        # use at own risk
# 
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
# 	package="AnnotationHub"
# 	moveFiles(package)

## ----eval=FALSE-----------------------------------------------------------------------------------
# 
#        path.expand(rappdirs::user_cache_dir(appname="AnnotationHub"))

## ----eval=FALSE-----------------------------------------------------------------------------------
# 
#     library(AnnotationHub)
#     package = "AnnotationHub"
# 
#     oldcache = path.expand(rappdirs::user_cache_dir(appname=package))
#     setAnnotationHubOption("CACHE", oldcache)
#     ah = AnnotationHub(localHub=TRUE)
#     ## removes old location and all resources
#     removeCache(ah, ask=FALSE)
# 
#     ## create the new default caching location
#     newcache = tools::R_user_dir(package, which="cache")
#     setAnnotationHubOption("CACHE", newcache)
#     ah = AnnotationHub()

