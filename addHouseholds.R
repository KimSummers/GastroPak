# addHouseholds
#
# Add household data to the plate counts
#
# file addHouseholds
#
# inputs
# 	countData     - Count data to have household data added
#   metaData      - Household data to add

# Version    Author       Date      Affiliation
# 1.00       J K Summers  09/01/24  Wellington Lab - School of Life Sciences - University of Warwick

addHouseholds <- function(countData, metaData) {

  for (iRow in 1:nrow(countData))
  {
    household <- paste("Household", metaData$Household[metaData$`Sample ID` ==
                                                         countData$SampleID[iRow]])

    if (iRow == 1)
    {
      households <- household
    }else
    {
      households <- rbind(households, household)
    }

  }

  countData <- cbind(countData, households)
  colnames(countData)[ncol(countData)] <- "Household"

  return(countData)
}
