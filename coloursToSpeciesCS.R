# coloursToSpeciesCS
#
# Add columns for species coubts based on colour counts on ChromoSelect agar
#
# file coloursToSpeciesCS
#
# inputs
# 	countData  - Initial Data with colour counts

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

coloursToSpeciesCS <- function(countData) {

  for (iRow in 1:nrow(countData))
  {
    eColiCount <- sum(countData$`Turquoise/Light Blue`[iRow],
                      countData$`Turquoise`[iRow], countData$`Violet/Purple`[iRow],
                      countData$`Purple/Voilet`[iRow], countData$`Dark Blue`[iRow])

    if (iRow == 1)
    {
      eColiCounts <- eColiCount
    }else
    {
      eColiCounts <- rbind(eColiCounts, eColiCount)
    }

  }

  countData <- cbind(countData, eColiCounts)
  rownames(countData) <- 1:nrow(countData)

  for (iRow in 1:nrow(countData))
  {
    coliformCount <- sum(countData$`Dark Pink`[iRow], countData$`Salmon to Red`[iRow])

    if (iRow == 1)
    {
      coliformCounts <- coliformCount
    }else
    {
      coliformCounts <- rbind(coliformCounts, coliformCount)
    }

  }

  countData <- cbind(countData, coliformCounts)

  return(countData)
}
