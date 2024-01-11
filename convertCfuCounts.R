# convertCfuCounts
#
# Convert counts in cfu form to raw data
#
# file convertCfuCounts
#
# inputs
# 	countData     - Initial Data to have counts averaged
#   numDilutions  - Number of replicates of data
#   convertCols   - The columns with counts to be converted

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

convertCfuCounts <- function(countData, numDilutions, convertCols) {

  for (iRow in 1:nrow(countData))
  {

    for (iDilution in 1:numDilutions)
    {

      if ((iRow %% numDilutions) == iDilution)
      {
        countData[iRow, convertCols] <- countData[iRow, convertCols] / 10^numDilutions
      }

    }

  }

  return(countData)
}
