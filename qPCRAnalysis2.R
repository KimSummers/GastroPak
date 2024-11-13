# qPCRAnalysis
#
# Read in qPCR data filter out sample ID, target, meanCount (gene copies) and SD
# where Ct for 2 or more values is undetermined or higher than the LOD meanCount
# is ND (not detected) and SD is NaN, for samples which are not ND, but 2 or more
# samples are undetermined or above the LOQ the mean count is D (detected) and
# the SD is again NaN
#
# file qPCRAnalysis.R
#
# inputs
# 	qPCRFileName          - file containing qPCRData data
#   dilutionFileName      - file containing dilution data
#   weightFile            - file containing weight data
# 	qPCRSumData           - existing summary data to add to
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  26/03/24  Wellington Lab - School of Life Sciences - University of Warwick
qPCRAnalysis2 <- function(qPCRFileName, qPCRSumData) {

  library(purrr)
  library(readxl)

  qPCRData <- read_excel(qPCRFileName, sheet = "Results")
  startRow <- which(qPCRData$`Block Type` == "Well")

  colnames(qPCRData) <- qPCRData[startRow, ]
  qPCRData <- as.data.frame(qPCRData)

  qPCRData <- qPCRData[(startRow + 1):nrow(qPCRData), ]

  if (is_empty(which(colnames(qPCRData) == "CT")))
  {
    colnames(qPCRData)[which(colnames(qPCRData) == "Cт")] <- "CT"
  }

  qPCRData$CT <- as.numeric(qPCRData$CT)

  if (is_empty(which(colnames(qPCRData) == "Well Position")))
  {
    qPCRSampleData <- qPCRData[which(qPCRData$Well == "C1"):nrow(qPCRData), ]
  }else
  {
    qPCRSampleData <- qPCRData[which(qPCRData$`Well Position` == "C1"):nrow(qPCRData), ]
  }

  qPCRSampleData <- qPCRSampleData[!is.na(qPCRSampleData$`Sample Name`), ]

  LOD <- as.numeric(qPCRData$LOD[1])
  LOQ <- as.numeric(qPCRData$LOQ[1])

  sampleNames <- unique(qPCRSampleData$`Sample Name`)

  for (iSample in 1:length(sampleNames))
  {
    sampleData <- qPCRSampleData[qPCRSampleData$`Sample Name` == sampleNames[iSample], ]

    if ((nrow(sampleData[is.na(sampleData$CT) | sampleData$CT > LOD, ])) > 1)
    {
      meanGeneCopy <- "ND"
      sdGeneCopy <- NaN
    }else
    {

      if (nrow(sampleData[sampleData$CT < LOQ, ]) < 2)
      {
        meanGeneCopy <- "D"
        sdGeneCopy <- NaN
      }else
      {
        meanGeneCopy <- as.numeric(sampleData$`Quantity Mean`[1])
        sdGeneCopy <- as.numeric(sampleData$`Quantity SD`[1])
      }

    }

    sampleRow <- sampleData[1, ]
    volDNA <- 35
    weightSample <- 250

    if (is.numeric(meanGeneCopy))
    {
      gcPerExtract <- meanGeneCopy * 100 * 50 / 3 / volDNA
      gcPerGOrPerMl <- gcPerExtract / weightSample
    }else
    {
      gcPerExtract <- meanGeneCopy
      gcPerGOrPerMl <- meanGeneCopy
    }

    if (is_empty(volDNA))
    {
      volDNA <- "Dilution Data Missing"
      gcPerGOrPerMl <- "ND"
      gcPerExtract <- "ND"
    }

    if (is_empty(weightSample))
    {
      weightSample <- "Extraction weight data missing"
      gcPerGOrPerMl <- "ND"
    }

    dataRow <- c(sampleRow$`Sample Name`, sampleRow$`Target Name`, meanGeneCopy,
                 sdGeneCopy, volDNA, weightSample, gcPerExtract, gcPerGOrPerMl)

    if (!is.null(qPCRSumData))
    {
      qPCRSumData <- rbind(qPCRSumData, dataRow)
    }else
    {
      qPCRSumData <- dataRow
    }

  }

  qPCRSumData <- as.data.frame(qPCRSumData)
  colnames(qPCRSumData) <- c("SampleID", "Target", "MeanGeneCopies",
                             "SdGeneCopies", "volumeDNA", "ExtractWeight",
                             "gcPerExtract", "gcPerGOrPerMl")

  rownames(qPCRSumData) <- 1:nrow(qPCRSumData)
  return(qPCRSumData)
}
