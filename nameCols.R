# nameCols
#
# Corrects column names
#
# file nameCols
#
# inputs
# 	countData     - Count data to have column names corrected

# Version    Author       Date      Affiliation
# 1.00       J K Summers  08/01/24  Wellington Lab - School of Life Sciences - University of Warwick

nameCols <- function(countData) {

  colnames(countData)[which(colnames(countData) == "Sample-ID")] <- "SampleID"
  colnames(countData)[which(colnames(countData) == "Sample ID")] <- "SampleID"
  colnames(countData)[which(colnames(countData) == "Sample")] <- "SampleID"
  colnames(countData)[which(colnames(countData) == "Sampling Code")] <- "SamplingCode"
  colnames(countData)[which(colnames(countData) == "Sampling Site")] <- "SamplingSite"
  colnames(countData)[which(colnames(countData) == "Location")] <- "SamplingSite"
  colnames(countData)[which(colnames(countData) == "Sample Type")] <- "SampleType"
  colnames(countData)[which(colnames(countData) == "Type")] <- "SampleType"
  colnames(countData)[which(colnames(countData) == "Light Pink")] <- "Klebsiella"
  colnames(countData)[which(colnames(countData) == "coliformCounts")] <- "coliforms"
  colnames(countData)[which(colnames(countData) == "eColiCounts")] <- "E.coli"
  colnames(countData)[which(colnames(countData) == "salmonellaCounts")] <- "Salmonella"
  colnames(countData)[which(colnames(countData) == "shigellaCounts")] <- "Shigella"
  colnames(countData)[which(colnames(countData) == "Test Country")] <- "TestCountry"
  colnames(countData)[which(colnames(countData) == "Country")] <- "TestCountry"
  colnames(countData)[which(colnames(countData) == "Gene copies in sample")] <- "MeanCfu"

  return(countData)
}
