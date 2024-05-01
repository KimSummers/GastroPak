# coloursToSpeciesBGA
#
# Add columns for species coubts based on colour counts on BG agar
#
# file coloursToSpeciesBGA
#
# inputs
# 	countData  - Initial Data with colour counts

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

coloursToSpeciesBGA <- function(countData) {

  salmonellaCounts <- NULL

  for (iRow in 1:nrow(countData))
  {
    salmonellaCount <- sum(countData$`Red to Pink`[iRow],
                            countData$`Red to Pinkish White`[iRow], na.rm = TRUE)
    salmonellaCounts <- rbind(salmonellaCounts, salmonellaCount)
  }

  countData <- cbind(countData, salmonellaCounts)
  rownames(countData) <- 1:nrow(countData)

  eColiCounts <- NULL

  for (iRow in 1:nrow(countData))
  {
    eColiCount <- sum(countData$`Partial Green`[iRow], countData$Yellow[iRow],
                      countData$`Yellow/ Yellowish green`[iRow], na.rm = TRUE)
    eColiCounts <- rbind(eColiCounts, eColiCount)
  }

  countData <- cbind(countData, eColiCounts)

  return(countData)
}
