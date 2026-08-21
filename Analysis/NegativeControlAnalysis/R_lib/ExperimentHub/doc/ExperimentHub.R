## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----library, message=FALSE---------------------------------------------------
library(ExperimentHub)

## ----ExperimentHub------------------------------------------------------------
eh = ExperimentHub()

## ----show---------------------------------------------------------------------
eh

## ----dataprovider-------------------------------------------------------------
head(unique(eh$dataprovider))

## ----species------------------------------------------------------------------
head(unique(eh$species))

## ----rdataclass---------------------------------------------------------------
head(unique(eh$rdataclass))

## ----alpine-------------------------------------------------------------------
apData <- query(eh, "alpineData")
apData

## ----show2--------------------------------------------------------------------
apData$preparerclass
df <- mcols(apData)

## ----length-------------------------------------------------------------------
length(eh)

## ----subset-------------------------------------------------------------------
mm <- query(eh, "mus musculus")
mm

## ----BiocHubsShiny, eval=FALSE------------------------------------------------
# BiocHubsShiny::BiocHubsShiny()

## ----echo=FALSE,fig.cap="BiocHubsShiny query with species term 'musculus'",out.width="100%"----
knitr::include_graphics("../man/figures/mm10_BiocHubsShiny.png")

## ----dm2----------------------------------------------------------------------
apData
apData["EH166"]

## ----dm3----------------------------------------------------------------------
apData[["EH166"]]

## ----show-2-------------------------------------------------------------------
eh

## ----snapshot-----------------------------------------------------------------
snapshotDate(eh)

## ----possibleDates------------------------------------------------------------
pd <- possibleDates(eh)
pd

## ----setdate, eval=FALSE------------------------------------------------------
# snapshotDate(ah) <- pd[1]

## ----eval=FALSE---------------------------------------------------------------
# 
#        # make sure you have permissions on the cache/files
#        # use at own risk
# 
# 
#        moveFiles<-function(package){
# 	   olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
# 	   newdir <- tools::R_user_dir(package, which="cache")
# 	   dir.create(path=newdir, recursive=TRUE)
# 	   files <- list.files(olddir, full.names =TRUE)
# 	   moveres <- vapply(files,
# 	       FUN=function(fl){
# 		   filename = basename(fl)
# 		   newname = file.path(newdir, filename)
# 		   file.rename(fl, newname)
# 	       },
# 	       FUN.VALUE = logical(1))
# 	   if(all(moveres)) unlink(olddir, recursive=TRUE)
#        }
# 
# 
#        package="ExperimentHub"
#        moveFiles(package)
# 

## ----eval=FALSE---------------------------------------------------------------
# path.expand(rappdirs::user_cache_dir(appname="ExperimentHub"))

## ----eval=FALSE---------------------------------------------------------------
# library(ExperimentHub)
# oldcache = path.expand(rappdirs::user_cache_dir(appname="ExperimentHub"))
# setExperimentHubOption("CACHE", oldcache)
# eh = ExperimentHub(localHub=TRUE)
# 
# ## removes old location and all resources
# removeCache(eh, ask=FALSE)
# 
# ## create the new default caching location
# newcache = tools::R_user_dir("ExperimentHub", which="cache")
# setExperimentHubOption("CACHE", newcache)
# eh = ExperimentHub()

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

