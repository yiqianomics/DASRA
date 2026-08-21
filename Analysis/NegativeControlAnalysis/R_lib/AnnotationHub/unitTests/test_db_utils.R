test_db_connections <- function(){

   hub <-  AnnotationHub()
   path <- hub@.db_path

   # test dbfile
   checkIdentical(path, dbfile(hub))

   conn <- dbconn(AnnotationHub())
   conn2 <-  AnnotationDbi::dbFileConnect(path)

   # test dbconn
   checkTrue(RSQLite::dbIsValid(conn))
   checkTrue(RSQLite::dbIsValid(conn2))
   checkIdentical(DBI::dbListTables(conn), DBI::dbListTables(conn2))

   DBI::dbDisconnect(conn)
   AnnotationHub:::.db_close(conn2)

   # test .db_uid0
   dates <- possibleDates(hub)
   ids1 <- AnnotationHub:::.db_uid0(path, dates[1])
   ids2 <- AnnotationHub:::.db_uid0(path, dates[length(dates)])
   # assumes more resources added than removed
   checkTrue(length(ids2) > length(ids1))

   # test .db_id
   checkIdentical(ids2,  AnnotationHub:::.db_uid(hub))

   # test index file
   indexfile <- AnnotationHub:::.db_index(hub)
   time1 <- file.mtime(indexfile)
   # this shouldn't recreate but just return existing
   indexfile2 <- AnnotationHub:::.db_create_index(hub)
   checkIdentical(indexfile, indexfile2)
   checkIdentical(time1, file.mtime(indexfile2))
   # same should just return
   indexfile3 <- AnnotationHub:::.db_index_file(hub)
   checkIdentical(unname(indexfile), indexfile3)
   checkIdentical(time1, file.mtime(indexfile3))
   index1 <- AnnotationHub:::.db_index_load(hub)
   index2 <- readRDS(indexfile)
   checkIdentical(index1, index2)

}


test_db_corrupt_and_fix <- function(){

    ## create hub to corrupt with multiple index files
    testhub = AnnotationHub(cache=tempfile(), ask=FALSE)
    cache = hubCache(testhub)
    bfc = BiocFileCache(cache)
    index_name = "annotationhub.index.rds"
    bfcnew(bfc, rname = index_name,  ext="_hub_index.rds")

    error_occured <- FALSE
    testhub2 <- withCallingHandlers({
        AnnotationHub(cache=cache)
    }, warning = function(cond) {
        error_occured <<- TRUE
        invokeRestart("muffleWarning")
    })
    checkTrue(error_occured)

    ## check basic information on original and 'new' hub are identical
    checkIdentical(hubCache(testhub), hubCache(testhub2))
    checkIdentical(AnnotationHub:::.db_uid(testhub), AnnotationHub:::.db_uid(testhub2))
    ## the index file was regenerated so it should not be the same
    checkTrue(!identical(unname((AnnotationHub:::.db_index(testhub))),
                         unname((AnnotationHub:::.db_index(testhub2)))))

    testhub3 = AnnotationHub(cache=cache)
    ## this should not regenerate and be the same
    checkTrue(identical(unname((AnnotationHub:::.db_index(testhub3))),
                        unname((AnnotationHub:::.db_index(testhub2)))))
    checkIdentical(hubCache(testhub3), hubCache(testhub2))
    checkIdentical(AnnotationHub:::.db_uid(testhub3), AnnotationHub:::.db_uid(testhub2))

}
