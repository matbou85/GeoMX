################################################################################
## This file does:                                                            ##
## 1. GeoMx codes                                                             ##
################################################################################

#----[ load packages ]------#
{
  library(NanoStringNCTools)
  library(GeomxTools)
  library(GeoMxWorkflows)
  library(knitr)
  library(dplyr)
  library(ggforce)
  library(networkD3)
  library(scales) 
  library(reshape2)
  library(cowplot)
  library(umap)
  library(Rtsne)
  library(pheatmap)
  library(ggrepel)
  library(ggplot2)
  library(tidyr)
  library(reshape2)
  library(Seurat)
  library(SpatialDecon)
  library(patchwork)
  library(standR)
  library(SpatialExperiment)
  library(ggalluvial)
  library(data.table)
  library(edgeR)
  library(limma)
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(vegan)
  
  # if(packageVersion("GeomxTools") < "2.1" & 
  #    packageVersion("GeoMxWorkflows") >= "1.0.1"){
  #   stop("GeomxTools and Workflow versions do not match. Please use the same version. 
  #     This workflow is meant to be used with most current version of packages. 
  #     If you are using an older version of Bioconductor please reinstall GeoMxWorkflows and use vignette(GeoMxWorkflows) instead")
  # }
  
  # if(packageVersion("GeomxTools") > "2.1" & 
  #    packageVersion("GeoMxWorkflows") <= "1.0.1"){
  #   stop("GeomxTools and Workflow versions do not match. 
  #          Please use the same version, see install instructions above.")
  
}

setwd("") ## Set wd accordingly

#----[ load data ]------#
{
  ## First re-defined columns without space in column names ##
  meta <- read_excel("GeoMx Worksheet_updated.xlsx")
  colnames(meta)
  setnames(meta, "Anatomical location", "Anatomical_location")
  setnames(meta, "Roi type", "Type")
  setnames(meta, "Mouse ID", "Mouse_ID")
  meta$Anatomical_location <- gsub("/", "_", meta$Anatomical_location)
  meta$Condition <- paste(meta$Anatomical_location, meta$Segment, meta$Type)
  meta$Anatomical_location <- gsub(" ", "_", meta$Anatomical_location)
  writexl::write_xlsx(meta, "GeoMx_meta.xlsx")
  
  DCCFiles <- dir(file.path("./DCC-20250328"), 
                  pattern = ".dcc$",
                  full.names = TRUE, recursive = TRUE)
  PKCFiles <- unzip(zipfile = dir(file.path("./pkcs"), 
                                  pattern = ".zip$",
                                  full.names = TRUE, recursive = TRUE))
  SampleAnnotationFile <-
    dir(file.path(""), 
        pattern = "GeoMx_meta.xlsx",
        full.names = TRUE, recursive = TRUE)
  
  # load data
  obj <-
    readNanoStringGeoMxSet(dccFiles = DCCFiles,
                           pkcFiles = PKCFiles,
                           phenoDataFile = SampleAnnotationFile,
                           phenoDataSheet = "Sheet1",
                           phenoDataDccColName = "Sample_ID",
                           protocolDataColNames = c("aoi", "roi"),
                           experimentDataColNames = c("panel"))
}


#----[ Modules Used ]------#
{
  pkcs <- annotation(obj) ## Set annotation information
  modules <- gsub(".pkc", "", pkcs)
  kable(data.frame(PKCs = pkcs, modules = modules))
}


#----obj#----[ Experimental design visualization  to look at the different types of samples and ROI/AOI segments ]------#
{
  sankeyCols <- c("source", "target", "value")
  
  link1 <- dplyr::count(pData(obj), `Type`, Injection)
  link2 <- dplyr::count(pData(obj),  Injection, `Anatomical_location`)
  link3 <- dplyr::count(pData(obj),  `Anatomical_location`, Segment)
  
  colnames(link1) <- sankeyCols
  colnames(link2) <- sankeyCols
  colnames(link3) <- sankeyCols
  
  links <- rbind(link1,link2,link3)
  nodes <- unique(data.frame(name=c(links$source, links$target)))
  
  # sankeyNetwork is 0 based, not 1 based
  links$source <- as.integer(match(links$source,nodes$name)-1)
  links$target <- as.integer(match(links$target,nodes$name)-1)
  
  
  sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                Target = "target", Value = "value", NodeID = "name",
                units = "TWh", fontSize = 12, nodeWidth = 30)
}



#----[ Replicate QC graphs from QC report ]------#
{
  # df2 <- exprs(obj)
  # Shift 0 counts to 1
  obj <- shiftCountsOne(obj, useDALogic = TRUE)
  
  # Default QC cutoffs are commented in () adjacent to the respective parameters
  # study-specific values were selected after visualizing the QC results in more
  # detail below
  QC_params <-
    list(minSegmentReads = 1000, # Minimum number of reads
         percentTrimmed = 80,    # Minimum % of reads trimmed
         percentStitched = 80,   # Minimum % of reads stitched
         percentAligned = 75,    # Minimum % of reads aligned
         percentSaturation = 50, # Minimum sequencing saturation
         minNegativeCount = 1,   # Minimum negative control counts
         maxNTCCount = 9000,     # Maximum counts observed in NTC well
         minNuclei = 20,         # Minimum # of nuclei estimated
         minArea = 1000)         # Minimum segment area
  obj <-
    setSegmentQCFlags(obj, 
                      qcCutoffs = QC_params)        
  
  # Collate QC Results
  protocolData(obj)[["QCFlags"]][[6]] <- NULL # I removed this column to not explicitely not remove ROIs with low counts for negative control probes
  QCResults <- protocolData(obj)[["QCFlags"]]
  flag_columns <- colnames(QCResults)
  
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns], na.rm = TRUE),
                           Warning = colSums(QCResults[, flag_columns], na.rm = TRUE))
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS", na.rm = TRUE),
      sum(QCResults[, "QCStatus"] == "WARNING", na.rm = TRUE))
  
  # Graphical summaries of QC statistics plot function
  QC_histogram <- function(assay_data = NULL,
                           annotation = NULL,
                           fill_by = NULL,
                           thr = NULL,
                           scale_trans = NULL) {
    plt <- ggplot(assay_data,
                  aes_string(x = paste0("unlist(`", annotation, "`)"),
                             fill = fill_by)) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept = thr, lty = "dashed", color = "black") +
      theme_bw() + guides(fill = "none") +
      facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
      labs(x = annotation, y = "Segments, #", title = annotation)
    if(!is.null(scale_trans)) {
      plt <- plt +
        scale_x_continuous(trans = scale_trans)
    }
    plt
  }

  df <- sData(obj)
  
  col_by <- "Segment"
  QC_histogram(sData(obj), "Trimmed (%)", col_by, 80)
  QC_histogram(sData(obj), "Raw", col_by, 80)
  QC_histogram(sData(obj), "Stitched (%)", col_by, 80)
  QC_histogram(sData(obj), "Aligned (%)", col_by, 75)
  QC_histogram(sData(obj), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")
  QC_histogram(sData(obj), "Area", col_by, 1000, scale_trans = "log10")
  QC_histogram(sData(obj), "Nuclei", col_by, 20)
  
  col_by <- "Type"
  QC_histogram(sData(obj), "Trimmed (%)", col_by, 80)
  QC_histogram(sData(obj), "Raw", col_by, 80)
  QC_histogram(sData(obj), "Stitched (%)", col_by, 80)
  QC_histogram(sData(obj), "Aligned (%)", col_by, 75)
  QC_histogram(sData(obj), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")
  QC_histogram(sData(obj), "Area", col_by, 1000, scale_trans = "log10")
  QC_histogram(sData(obj), "Nuclei", col_by, 20)
  
  col_by <- "Injection"
  QC_histogram(sData(obj), "Trimmed (%)", col_by, 80)
  QC_histogram(sData(obj), "Raw", col_by, 80)
  QC_histogram(sData(obj), "Stitched (%)", col_by, 80)
  QC_histogram(sData(obj), "Aligned (%)", col_by, 75)
  QC_histogram(sData(obj), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")
  QC_histogram(sData(obj), "Area", col_by, 1000, scale_trans = "log10")
  QC_histogram(sData(obj), "Nuclei", col_by, 20)
  
  col_by <- "Anatomical_location"
  QC_histogram(sData(obj), "Trimmed (%)", col_by, 80)
  QC_histogram(sData(obj), "Raw", col_by, 80)
  QC_histogram(sData(obj), "Stitched (%)", col_by, 80)
  QC_histogram(sData(obj), "Aligned (%)", col_by, 75)
  QC_histogram(sData(obj), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")
  QC_histogram(sData(obj), "Area", col_by, 1000, scale_trans = "log10")
  QC_histogram(sData(obj), "Nuclei", col_by, 20)

  # calculate the negative geometric means for each module
  negativeGeoMeans <- 
    esBy(negativeControlSubset(obj), 
         GROUP = "Module", 
         FUN = function(x) { 
           assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
  protocolData(obj)[["NegGeoMean"]] <- negativeGeoMeans
  
  # explicitly copy the Negative geoMeans from sData to pData
  negCols <- paste0("NegGeoMean_", modules)
  pData(obj)[, negCols] <- sData(obj)[["NegGeoMean"]]
  for(ann in negCols) {
    plt <- QC_histogram(pData(obj), ann, col_by, 2, scale_trans = "log10")
    print(plt)
  }
  
  # detatch neg_geomean columns ahead of aggregateCounts call
  pData(obj) <- pData(obj)[, !colnames(pData(obj)) %in% negCols]
  
  # show all NTC values, Freq = # of Segments with a given NTC count:
  kable(table(NTC_Count = sData(obj)$NTC),
        col.names = c("NTC Count", "# of Segments"))
  
  kable(QC_Summary, caption = "QC Summary Table for each Segment")
  
  ## Exclude warnings
  obj <- obj[, QCResults$QCStatus == "PASS"]
  
  # Subsetting our dataset has removed samples which did not pass QC
  dim(obj)
  #> Features  Samples 
  #>    20465      88
  
}



#----[ Probe QC and remove outliers ]------#
{
  ## Set Probe QC Flags
  # Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
  # FALSE if you do not want to remove local outliers
  obj <- setBioProbeQCFlags(obj, 
                                 qcCutoffs = list(minProbeRatio = 0.1,
                                                  percentFailGrubbs = 20), 
                                 removeLocalOutliers = TRUE)
  
  ProbeQCResults <- fData(obj)[["QCFlags"]]
  
  # Define QC table for Probe QC
  qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                      Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                      Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                  & !ProbeQCResults$GlobalGrubbsOutlier))
  
  ## Exclude Outlier Probes
  #Subset object to exclude all that did not pass Ratio & Global testing
  ProbeQCPassed <- 
    subset(obj, 
           fData(obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
             fData(obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
  dim(ProbeQCPassed)
  #> Features  Samples 
  #>    20465      88
  obj <- ProbeQCPassed 
}



#----[ Create Gene-level Count Data ]------#
{
  # Check how many unique targets the object has
  length(unique(featureData(obj)[["TargetName"]]))
  #> [1] 20256
  
  # collapse to targets
  target_obj <- aggregateCounts(obj)
  dim(target_obj)
  #> Features  Samples 
  #>    20256      88
  exprs(target_obj)[1:5, 1:2]
  #> DSP-1001660039737-B-A02.dcc DSP-1001660039737-B-A03.dcc
  #> Trappc2b                                1                           3
  #> Sanbr                                   1                           1
  #> 0610010K14Rik                           1                           2
  #> Ncbp2as2                                1                           4
  #> 0610030E20Rik                           1                           1
}





#----[ Limit of Quantification (LOQ) ]------#
{
  # Define LOQ SD threshold and minimum value
  cutoff <- 2
  minLOQ <- 2
  
  # Calculate LOQ per module tested
  LOQ <- data.frame(row.names = colnames(target_obj))
  for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_obj)))) {
      LOQ[, module] <-
        pmax(minLOQ,
             pData(target_obj)[, vars[1]] * 
               pData(target_obj)[, vars[2]] ^ cutoff)
    }
  }
  pData(target_obj)$LOQ <- LOQ
}

# df <- pData(target_obj)
# write.table(as.data.frame(df), "metadata.txt", sep = "\t")

#----[ Filtering LOQ ]------#
{
  LOQ_Mat <- c()
  for(module in modules) {
    ind <- fData(target_obj)$Module == module
    Mat_i <- t(esApply(target_obj[ind, ], MARGIN = 1,
                       FUN = function(x) {
                         x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
  }
  # ensure ordering since this is stored outside of the geomxSet
  LOQ_Mat <- LOQ_Mat[fData(target_obj)$TargetName, ]
  dim(target_obj)
  
}



#----[ Filtering Segment ]------#
{
  # Save detection rate information to pheno data
  pData(target_obj)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
  pData(target_obj)$GeneDetectionRate <-
    pData(target_obj)$GenesDetected / nrow(target_obj)
  
  # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
  pData(target_obj)$DetectionThreshold <- 
    cut(pData(target_obj)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
  
  # stacked bar plot of different cut points (1%, 5%, 10%, 15%)
  ggplot(pData(target_obj),
         aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = Segment)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")
  #> Warning: The dot-dot notation (`..count..`) was deprecated in ggplot2 3.4.0.
  #> ℹ Please use `after_stat(count)` instead.
  #> This warning is displayed once every 8 hours.
  #> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
  #> generated.
  
  # cut percent genes detected at 1, 5, 10, 15
  kable(table(pData(target_obj)$DetectionThreshold,
              pData(target_obj)$Segment))
  
  target_obj <-
    target_obj[, pData(target_obj)$GeneDetectionRate >= .001]
  
  dim(target_obj)
  #> Features  Samples 
  #>    20256      83

}

# df3 <- exprs(target_obj)

#----[ Filtering Gene detection ]------#
{
  # Calculate detection rate:
  LOQ_Mat <- LOQ_Mat[, colnames(target_obj)]
  fData(target_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
  fData(target_obj)$DetectionRate <-
    fData(target_obj)$DetectedSegments / nrow(pData(target_obj))
  
  # Gene of interest detection table
  temp <- fData(target_obj)
  goi <- c("Col1a1", "Mmp9", "Mmp14", "Cd68", "Ctsk")
  
  goi_df <- data.frame(
    Gene = goi,
    Number = fData(target_obj)[goi, "DetectedSegments"],
    DetectionRate = percent(fData(target_obj)[goi, "DetectionRate"]))
  goi_df
  # Gene Number DetectionRate
  # 1 Col1a1     58         71.6%
  # 2   Mmp9     54         66.7%
  # 3  Mmp14     43         53.1%
  # 4   Cd68     29         35.8%
}


#----[ Filtering Genes ]------#
{
  # Plot detection rate:
  plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
  plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_obj)$DetectionRate >= x)}))
  plot_detect$Rate <- plot_detect$Number / nrow(fData(target_obj))
  rownames(plot_detect) <- plot_detect$Freq
  
  ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")
  
  # Subset to target genes detected in at least 10% of the samples.
  #   Also manually include the negative control probe, for downstream use
  negativeProbefData <- subset(fData(target_obj), CodeClass == "Negative")
  neg_probes <- unique(negativeProbefData$TargetName)
  target_obj <- 
    target_obj[fData(target_obj)$DetectionRate >= 0.1 |
                      fData(target_obj)$TargetName %in% neg_probes, ]
  dim(target_obj)
  #> Features  Samples 
  #>    2896      83
  
  # retain only detected genes of interest
  goi <- goi[goi %in% rownames(target_obj)]
}

# df <- exprs(target_obj)


#----[ Normalization ]------#
{
  # Graph Q3 value vs negGeoMean of Negatives
  ann_of_interest <- "Segment"
  Stat_data <- 
    data.frame(row.names = colnames(exprs(target_obj)),
               Segment = colnames(exprs(target_obj)),
               Annotation = pData(target_obj)[, ann_of_interest],
               Q3 = unlist(apply(exprs(target_obj), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = exprs(target_obj)[neg_probes, ])
  Stat_data_m <- reshape2::melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  plt1 <- ggplot(Stat_data_m,
                 aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")
  
  plt2 <- ggplot(Stat_data,
                 aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
  
  plt3 <- ggplot(Stat_data,
                 aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
  
  btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                       rel_widths = c(0.43,0.57))
  plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
  
  # Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
  target_obj <- normalize(target_obj ,
                               norm_method = "quant", 
                               desiredQuantile = .75,
                               toElt = "q_norm")
  
  # Background normalization for WTA/CTA without custom spike-in
  target_obj <- normalize(target_obj ,
                               norm_method = "neg", 
                               fromElt = "exprs",
                               toElt = "neg_norm")
  
  # visualize the first 10 segments with each normalization method
  boxplot(exprs(target_obj)[,1:10],
          col = "#9EDAE5", main = "Raw Counts",
          log = "y", names = 1:10, xlab = "Segment",
          ylab = "Counts, Raw")
  
  boxplot(assayDataElement(target_obj[,1:10], elt = "q_norm"),
          col = "#2CA02C", main = "Q3 Norm Counts",
          log = "y", names = 1:10, xlab = "Segment",
          ylab = "Counts, Q3 Normalized")
  
  boxplot(assayDataElement(target_obj[,1:10], elt = "neg_norm"),
          col = "#FF7F0E", main = "Neg Norm Counts",
          log = "y", names = 1:10, xlab = "Segment",
          ylab = "Counts, Neg. Normalized")
}


# objSPE <- as.SpatialExperiment(target_obj, normData = "neg_norm")
# data.frame(head(colData(objSPE)))
# assayNames(objSPE)



#----[ Transform object to be read for standR analysis ]------#
{
  setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
  counts <- target_obj@assayData[["exprs"]] |> as.data.frame()
  counts$TargetName <- rownames(counts)
  # head(counts)[,1:5]
  metadata <- pData(target_obj) |> as.data.frame()
  metadata$SegmentDisplayName <- rownames(metadata)
  # head(metadata)
  genemeta <- target_obj@featureData@data |> as.data.frame()
  # head(genemeta)[,1:5]
  all(genemeta$TargetName == counts$TargetName)
  objSPE <- standR::readGeoMx(counts, metadata, genemeta, 
                      coord.colnames = c("ROI Coordinate X", "ROI Coordinate Y"))
  assayNames(objSPE)
  assay(objSPE, "counts")[1:5,1:5]
  assay(objSPE, "logcounts")[1:5,1:5]
  colData(objSPE)[1:5,1:5]
  rowData(objSPE)[1:5,1:5]
  metadata(objSPE)$NegProbes[,1:5]
  colData(objSPE)$QCFlags
  
  ## set 2 counts as a threshold  in 90% of the samples
  objSPE <- addPerROIQC(objSPE, min_count = 2, sample_fraction = 0.9, rm_genes = TRUE)
  dim(objSPE)
  metadata(objSPE) |> names()  
  plotGeneQC(objSPE, ordannots = "Segment", col = Segment, point_size = 2)
  
  ## set 3 counts as a threshold  in 90% of the samples
  objSPE <- addPerROIQC(objSPE, min_count = 3, sample_fraction = 0.9, rm_genes = TRUE)
  dim(objSPE)
  metadata(objSPE) |> names()
  
  ## The percentage of lowly expressed genes per sample is rather high
  
  ## Set threshold as 10 nuclei count
  plotROIQC(objSPE, x_threshold = 1, color = 'Slide name',   x_axis = "Nuclei") 
  ## The number of nuclei count is also low
  qc <- colData(objSPE)$Nuclei > 1
  table(qc)
  
  objSPE <- objSPE[, qc]
  
  plotROIQC(objSPE,  x_axis = "Area", x_lab = "AreaSize", 
            y_axis = "lib_size", y_lab = "Library size", 
            col = 'Slide name')
  
  ## Plot raw counts
  plotRLExpr(objSPE) 
  ## Plot normalised counts grouped by segments
  plotRLExpr(objSPE, ordannots = "Segment", assay = 2, color = Segment) 
  ## Plot normalised counts grouped by segments
  plotRLExpr(objSPE, ordannots = "Anatomical_location", assay = 2, 
             color = Anatomical_location) 
  ## Plot normalised counts grouped by segments
  plotRLExpr(objSPE, ordannots = "Type", assay = 2, color = Type) 
  
  ## PCA
  set.seed(1985)
  objSPE <- scater::runPCA(objSPE)
  pca_results <- reducedDim(objSPE, "PCA")
  plotScreePCA(objSPE, precomputed = pca_results)
  plotPairPCA(objSPE, col = Type, precomputed = pca_results, n_dimension = 4)
  plotPairPCA(objSPE, col = Segment, precomputed = pca_results, n_dimension = 4)
  plotPairPCA(objSPE, col = Injection, precomputed = pca_results, n_dimension = 4)
  drawPCA(objSPE, precomputed = pca_results, dims = c(1, 3), col = Type)
  drawPCA(objSPE, precomputed = pca_results, col = Segment)
  ## Top 10 genes influencing PCA
  plotPCAbiplot(objSPE, n_loadings = 10, dims = c(1, 3), 
                precomputed = pca_results, col = Type) 
  
  ## MDS
  standR::plotMDS(objSPE, assay = 2, color = Type)
  
  ## UMAP
  set.seed(1985)
  objSPE <- scater::runUMAP(objSPE, dimred = "PCA")
  plotDR(objSPE, dimred = "UMAP", col = Mouse_ID)
  plotDR(objSPE, dimred = "UMAP", col = Type)
  plotDR(objSPE, dimred = "UMAP", col = Injection) 
  plotDR(objSPE, dimred = "UMAP", col = Segment) 
  
  plotDR(objSPE, dimred = "UMAP", col = Injection) +
    geom_text_repel(
      aes(label = rownames(colData(objSPE))),
      max.overlaps = Inf,
      size=2
    )
  
  ## The TMM normalization seems to be a better option for these data
  objSPE_tmm <- geomxNorm(objSPE, method = "TMM")
  plotRLExpr(objSPE_tmm, assay = 2, color = Type) + ggtitle("TMM")

  saveRDS(objSPE_tmm, file = "GeoMx_SpatialExperiment.RDS")
  
}



#----[ Remove unwanted batch effects ]------#
{
  colData(objSPE)
  objSPE <- findNCGs(objSPE, batch_name = "slide name", top_n = 300)
  metadata(objSPE) |> names()
  
  for(i in seq(5)){
    spe_ruv <- geomxBatchCorrection(objSPE, factors = "Type", 
                                    NCGs = metadata(objSPE)$NCGs, k = i)
    
    print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Type, title = paste0("k = ", i)))
    
  }
  
  spe_ruv <- geomxBatchCorrection(objSPE, factors = "Type", 
                                  NCGs = metadata(objSPE)$NCGs, k = 3)
  set.seed(1985)
  spe_ruv <- scater::runPCA(spe_ruv)
  pca_results_ruv <- reducedDim(spe_ruv, "PCA")
  plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = Type, title = "RUV4, k = 3", n_dimension = 4)
  
  tmp <- data.frame()
  tmp <- cbind(pca_results_ruv[,1], pca_results_ruv[,2])
  group <- colData(spe_ruv) |> as.data.frame()
  group <- group$Type
  df <- cbind(tmp, group) |> as.data.frame()
  test_stats <- adonis2(tmp ~ group, method = "eu")
  drawPCA(spe_ruv, precomputed = pca_results_ruv, col = Type) %>%
    + annotate(geom = "text", 
               x = -Inf, y = Inf, 
               hjust = 0, vjust = 1.5,
               label = paste0("PERMANOVA test: F = ", 
                              format(signif(test_stats$F[[1]], digits = 2), nsmall = 2), ", P = ", 
                              format(signif(test_stats$`Pr(>F)`[[1]], digits = 2), nsmall = 2)),
               color = "black", 
               size = 6)
  
  
  plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = SlideName, title = "RUV4, k = 2", n_dimension = 4)
  
}













#----[ Unsupervised Analysis ]------#
{
  # update defaults for umap to contain a stable random_state (seed)
  custom_umap <- umap::umap.defaults
  custom_umap$random_state <- 42
  # run UMAP
  umap_out <-
    umap(t(log2(assayDataElement(target_obj , elt = "q_norm"))),  
         config = custom_umap)
  #> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
  #> Also defined by 'spam'
  pData(target_obj)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
  ggplot(pData(target_obj),
         aes(x = UMAP1, y = UMAP2, color = Segment, shape = Segment)) +
    geom_point(size = 3) +
    theme_bw()
  
  # run tSNE
  set.seed(1985) # set the seed for tSNE as well
  tsne_out <-
    Rtsne(t(log2(assayDataElement(target_obj , elt = "q_norm"))),
          perplexity = ncol(target_obj)*.15)
  pData(target_obj)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
  ggplot(pData(target_obj),
         aes(x = tSNE1, y = tSNE2, color = Segment, shape = Segment)) +
    geom_point(size = 3) +
    theme_bw()
}


#----[ Clustering high cv genes ]------#
{
  # create a log2 transform of the data for analysis
  assayDataElement(object = target_obj, elt = "log_q") <-
    assayDataApply(target_obj, 2, FUN = log, base = 2, elt = "q_norm")
  
  # create CV function
  calc_CV <- function(x) {sd(x) / mean(x)}
  CV_dat <- assayDataApply(target_obj,
                           elt = "log_q", MARGIN = 1, calc_CV)
  # show the highest CD genes and their CV values
  sort(CV_dat, decreasing = TRUE)[1:5]
  #>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
  #> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216
  
  # Identify genes in the top 3rd of the CV values
  GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
  pheatmap(assayDataElement(target_obj[GOI, ], elt = "log_q"),
           scale = "row", 
           show_rownames = FALSE, show_colnames = FALSE,
           border_color = NA,
           clustering_method = "average",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           breaks = seq(-3, 3, 0.05),
           color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
           annotation_col = 
             pData(target_obj)[, c("class", "segment", "region")])
}



#----[ Differential Expression Within Slide Analysis ]------#
{
  # convert test variables to factors
  pData(target_demoData)$testRegion <-
    factor(pData(target_demoData)$region, c("glomerulus", "tubule"))
  pData(target_demoData)[["slide"]] <-
    factor(pData(target_demoData)[["slide name"]])
  assayDataElement(object = target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")
  
  # run LMM:
  # formula follows conventions defined by the lme4 package
  results <- c()
  for(status in c("DKD", "normal")) {
    ind <- pData(target_demoData)$class == status
    mixedOutmc <-
      mixedModelDE(target_demoData[, ind],
                   elt = "log_q",
                   modelFormula = ~ testRegion + (1 + testRegion | slide),
                   groupVar = "testRegion",
                   nCores = parallel::detectCores(),
                   multiCore = FALSE)
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <-
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results <- rbind(results, r_test)
  }
}



#----[ Interpreting the results table ]------#
{
  kable(subset(results, Gene %in% goi & Subset == "normal"), digits = 3,
        caption = "DE results for Genes of Interest",
        align = "lc", row.names = FALSE)
}



#----[ Differential Expression between Slide Analysis ]------#
{
  # convert test variables to factors
  pData(target_demoData)$testClass <-
    factor(pData(target_demoData)$class, c("normal", "DKD"))
  
  # run LMM:
  # formula follows conventions defined by the lme4 package
  results2 <- c()
  for(region in c("glomerulus", "tubule")) {
    ind <- pData(target_demoData)$region == region
    mixedOutmc <-
      mixedModelDE(target_demoData[, ind],
                   elt = "log_q",
                   modelFormula = ~ testClass + (1 | slide),
                   groupVar = "testClass",
                   nCores = parallel::detectCores(),
                   multiCore = FALSE)
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <-
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- region
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results2 <- rbind(results2, r_test)
  }
  
  
  kable(subset(results2, Gene %in% goi & Subset == "tubule"), digits = 3,
        caption = "DE results for Genes of Interest",
        align = "lc", row.names = FALSE)
}



#----[ Visualizing DE Genes ]------#
{
  ## Volcano plots
  
  # Categorize Results based on P-value & FDR for plotting
  results$Color <- "NS or FC < 0.5"
  results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results$Color[results$FDR < 0.05] <- "FDR < 0.05"
  results$Color[results$FDR < 0.001] <- "FDR < 0.001"
  results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
  results$Color <- factor(results$Color,
                          levels = c("NS or FC < 0.5", "P < 0.05",
                                     "FDR < 0.05", "FDR < 0.001"))
  
  # pick top genes for either side of volcano to label
  # order genes for convenience:
  results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
  top_g <- c()
  for(cond in c("DKD", "normal")) {
    ind <- results$Subset == cond
    top_g <- c(top_g,
               results[ind, 'Gene'][
                 order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results[ind, 'Gene'][
                 order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
  }
  top_g <- unique(top_g)
  results <- results[, -1*ncol(results)] # remove invert_P from matrix
  
  # Graph results
  ggplot(results,
         aes(x = Estimate, y = -log10(`Pr(>|t|)`),
             color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = "Enriched in Tubules <- log2(FC) -> Enriched in Glomeruli",
         y = "Significance, -log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.001),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    facet_wrap(~Subset, scales = "free_y")
  
  ## Plotting Genes of Interest
  kable(subset(results, Gene %in% c('PDHA1','ITGB1')), row.names = FALSE)
  
  # show expression for a single target: PDHA1
  ggplot(pData(target_demoData),
         aes(x = region, fill = region,
             y = as.numeric(assayDataElement(target_demoData["PDHA1", ],
                                             elt = "q_norm")))) +
    geom_violin() +
    geom_jitter(width = .2) +
    labs(y = "PDHA1 Expression") +
    scale_y_continuous(trans = "log2") +
    facet_wrap(~class) +
    theme_bw()
  
  glom <- pData(target_demoData)$region == "glomerulus"
  
  # show expression of PDHA1 vs ITGB1
  ggplot(pData(target_demoData),
         aes(x = as.numeric(assayDataElement(target_demoData["PDHA1", ],
                                             elt = "q_norm")),
             y = as.numeric(assayDataElement(target_demoData["ITGB1", ],
                                             elt = "q_norm")),
             color = region)) +
    geom_vline(xintercept =
                 max(assayDataElement(target_demoData["PDHA1", glom],
                                      elt = "q_norm")),
               lty = "dashed", col = "darkgray") +
    geom_hline(yintercept =
                 max(assayDataElement(target_demoData["ITGB1", !glom],
                                      elt = "q_norm")),
               lty = "dashed", col = "darkgray") +
    geom_point(size = 3) +
    theme_bw() +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    labs(x = "PDHA1 Expression", y = "ITGB1 Expression") +
    facet_wrap(~class)
  
  ## Heatmap of Significant Genes
  # select top significant genes based on significance, plot with pheatmap
  GOI <- unique(subset(results, `FDR` < 0.001)$Gene)
  pheatmap(log2(assayDataElement(target_demoData[GOI, ], elt = "q_norm")),
           scale = "row",
           show_rownames = FALSE, show_colnames = FALSE,
           border_color = NA,
           clustering_method = "average",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           cutree_cols = 2, cutree_rows = 2,
           breaks = seq(-3, 3, 0.05),
           color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
           annotation_col = pData(target_demoData)[, c("region", "class")])
  
}



#----[ Visualizing DE Genes ]------#
{
  demoSeurat <- as.Seurat(target_demoData, normData = "exprs")
  demoSeurat
  head(demoSeurat, 3) # most important ROI metadata
  demoSeurat@misc[1:8] # experiment data
  head(demoSeurat@misc$sequencingMetrics) # sequencing metrics
  head(demoSeurat@misc$QCMetrics$QCFlags) # QC metrics
  head(demoSeurat@assays$GeoMx@meta.data) # gene metadata
  VlnPlot(demoSeurat, features = "nCount_GeoMx", pt.size = 0.1)
  demoSeurat <- as.Seurat(target_demoData, normData = "q_norm", ident = "region")
  VlnPlot(demoSeurat, features = "nCount_GeoMx", pt.size = 0.1)
  
  }


#----[ Normalization ]------#
{
  demoSeurat <- FindVariableFeatures(demoSeurat)
  demoSeurat <- ScaleData(demoSeurat)
  demoSeurat <- RunPCA(demoSeurat, assay = "GeoMx", verbose = FALSE)
  demoSeurat <- FindNeighbors(demoSeurat, reduction = "pca", dims = seq_len(30))
  demoSeurat <- FindClusters(demoSeurat, verbose = FALSE)
  demoSeurat <- RunUMAP(demoSeurat, reduction = "pca", dims = seq_len(30))
  
  DimPlot(demoSeurat, reduction = "umap", label = TRUE, group.by = "region")
}



#----[ In depth ]------#
{
  data("nsclc", package = "SpatialDecon")
  nsclc
  dim(nsclc)
  data.frame(exprs(nsclc)[seq_len(5), seq_len(5)])
  head(pData(nsclc))
  head(pData(target_demoData))
  
  nsclcSeurat <- as.Seurat(nsclc, normData = "q_norm", ident = "aoi", 
                           coordinates = c("x", "y"))
  
  nsclcSeurat
  
  
}


















