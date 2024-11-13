# fixMissingData
#
# Fix gaps in isolations data
#
# file fixMissingData
#
# inputs
# 	fixData  - Data to have missing values fixed
#   fixReps  - Number of replicates of data

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

fixMissingData <- function(fixData, fixReps) {

fixCount <- 3 * fixReps

for (iRow in (1:(nrow(fixData) / fixCount)))
{
  fixData$`Sampling Site`[((iRow - 1) * fixCount + 2):((iRow - 1) * fixCount + fixCount)] <-
    fixData$`Sampling Site`[((iRow - 1) * fixCount + 1)]
  fixData$`Sample Type`[((iRow - 1) * fixCount + 2):((iRow - 1) * fixCount + fixCount)] <-
    fixData$`Sample Type`[((iRow - 1) * fixCount + 1)]
}

for (iRow in (1:(nrow(fixData) / fixCount)))
{
  fixData$`SamplingSite`[((iRow - 1) * fixCount + 2):((iRow - 1) * fixCount + fixCount)] <-
    fixData$`SamplingSite`[((iRow - 1) * fixCount + 1)]
  fixData$`SampleType`[((iRow - 1) * fixCount + 2):((iRow - 1) * fixCount + fixCount)] <-
    fixData$`SampleType`[((iRow - 1) * fixCount + 1)]
}

for (iRow in (1:(nrow(fixData) / fixReps)))
{
  fixData$`SamplingCode`[((iRow - 1) * fixReps + 2):((iRow - 1) * fixReps + fixReps)] <-
    fixData$`SamplingCode`[((iRow - 1) * fixReps + 1)]
  fixData$Replicate[((iRow - 1) * fixReps + 2):((iRow - 1) * fixReps + fixReps)] <-
    fixData$Replicate[((iRow - 1) * fixReps + 1)]
  fixData$`SampleID`[((iRow - 1) * fixReps + 2):((iRow - 1) * fixReps + fixReps)] <-
    fixData$`SampleID`[((iRow - 1) * fixReps + 1)]
}

return(fixData)
}
