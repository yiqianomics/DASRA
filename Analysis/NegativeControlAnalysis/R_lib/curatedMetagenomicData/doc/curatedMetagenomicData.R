## ----include = FALSE----------------------------------------------------------
library(curatedMetagenomicData)

## ----eval = FALSE-------------------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("curatedMetagenomicData")

## ----message = FALSE----------------------------------------------------------
library(dplyr)
library(DT)

## ----collapse = TRUE----------------------------------------------------------
sampleMetadata |>
    filter(study_name == "AsnicarF_2017") |>
    select(where(~ !any(is.na(.x)))) |>
    slice(1:10) |>
    select(1:10) |>
    datatable(options = list(dom = "t"), extensions = "Responsive")

## ----collapse = TRUE----------------------------------------------------------
curatedMetagenomicData("AsnicarF_20.+")

## ----collapse = TRUE, message = FALSE-----------------------------------------
curatedMetagenomicData("AsnicarF_2017.relative_abundance", dryrun = FALSE, rownames = "short")

## ----collapse = TRUE, message = FALSE-----------------------------------------
curatedMetagenomicData("AsnicarF_20.+.relative_abundance", dryrun = FALSE, counts = TRUE, rownames = "short")

## ----collapse = TRUE, message = FALSE-----------------------------------------
curatedMetagenomicData("AsnicarF_20.+.marker_abundance", dryrun = FALSE) |>
    mergeData()

## ----collapse = TRUE, message = FALSE-----------------------------------------
curatedMetagenomicData("AsnicarF_20.+.pathway_abundance", dryrun = FALSE) |>
    mergeData()

## ----collapse = TRUE, message = FALSE-----------------------------------------
curatedMetagenomicData("AsnicarF_20.+.relative_abundance", dryrun = FALSE, rownames = "short") |>
    mergeData()

## ----collapse = TRUE, message = FALSE-----------------------------------------
sampleMetadata |>
    filter(age >= 18) |>
    filter(!is.na(alcohol)) |>
    filter(body_site == "stool") |>
    select(where(~ !all(is.na(.x)))) |>
    returnSamples("relative_abundance", rownames = "short")

## ----message = FALSE----------------------------------------------------------
library(stringr)
library(mia)
library(scater)
library(vegan)

## ----collapse = TRUE, message = FALSE-----------------------------------------
alcoholStudy <-
    filter(sampleMetadata, age >= 18) |>
    filter(!is.na(alcohol)) |>
    filter(body_site == "stool") |>
    select(where(~ !all(is.na(.x)))) |>
    returnSamples("relative_abundance", rownames = "short")

## ----collapse = TRUE, message = FALSE-----------------------------------------
colData(alcoholStudy) <-
    colData(alcoholStudy) |>
    as.data.frame() |>
    mutate(alcohol = str_replace_all(alcohol, "no", "No")) |>
    mutate(alcohol = str_replace_all(alcohol, "yes", "Yes")) |>
    DataFrame()

## ----collapse = TRUE, message = FALSE-----------------------------------------
altExps(alcoholStudy) <-
    splitByRanks(alcoholStudy)

## ----collapse = TRUE, fig.cap = "Alpha Diversity – Shannon Index (H')"--------
alcoholStudy |>
    estimateDiversity(assay.type = "relative_abundance", index = "shannon") |>
    plotColData(x = "alcohol", y = "shannon", colour_by = "alcohol", shape_by = "alcohol") +
    labs(x = "Alcohol", y = "Alpha Diversity (H')") +
    guides(colour = guide_legend(title = "Alcohol"), shape = guide_legend(title = "Alcohol")) +
    theme(legend.position = "none")

## ----collapse = TRUE, fig.cap = "Beta Diversity – Bray–Curtis PCoA"-----------
alcoholStudy |>
    runMDS(FUN = vegdist, method = "bray", exprs_values = "relative_abundance", altexp = "genus", name = "BrayCurtis") |>
    plotReducedDim("BrayCurtis", colour_by = "alcohol", shape_by = "alcohol") +
    labs(x = "PCo 1", y = "PCo 2") +
    guides(colour = guide_legend(title = "Alcohol"), shape = guide_legend(title = "Alcohol")) +
    theme(legend.position = c(0.90, 0.85))

## ----collapse = TRUE, fig.cap = "Beta Diversity – UMAP (Uniform Manifold Approximation and Projection)"----
alcoholStudy |>
    runUMAP(exprs_values = "relative_abundance", altexp = "genus", name = "UMAP") |>
    plotReducedDim("UMAP", colour_by = "alcohol", shape_by = "alcohol") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(colour = guide_legend(title = "Alcohol"), shape = guide_legend(title = "Alcohol")) +
    theme(legend.position = c(0.90, 0.85))

## ----eval = FALSE-------------------------------------------------------------
# convertToPhyloseq(alcoholStudy, assay.type = "relative_abundance")

## ----echo = FALSE-------------------------------------------------------------
utils::sessionInfo()

