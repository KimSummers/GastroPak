# coloursToSpeciesSSA
#
# Add columns for species coubts based on colour counts on SS agar
#
# file coloursToSpeciesSSA
#
# inputs
# 	countData  - Initial Data with colour counts

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

coloursToSpeciesSSA <- function(countData) {

  eColiCounts <- NULL

  for (iRow in 1:nrow(countData))
  {
    eColiCount <- sum(countData$`Pink, Bile precipitate`[iRow],
                      countData$`Pink Bile Precipitate`[iRow],
                      countData$`Red/Orange`[iRow],
                      countData$`Red to Pink`[iRow])

    eColiCounts <- rbind(eColiCounts, eColiCount)
  }

  countData <- cbind(countData, eColiCounts)
  rownames(countData) <- 1:nrow(countData)

  shigellaCounts <- NULL

  for (iRow in 1:nrow(countData))
  {
    shigellaCount <- sum(countData$Colorless[iRow], countData$Yellow[iRow])

    shigellaCounts <- rbind(shigellaCounts, shigellaCount)
  }

  countData <- cbind(countData, shigellaCounts)

  salmonellaCounts <- NULL

  for (iRow in 1:nrow(countData))
  {
    salmonellaCount <- sum(countData$`Colorless with Black center`[iRow],
                           countData$`Colorless with Brown center`[iRow],
                           countData$`Colorless with Yellow center`[iRow],
                           countData$`Colorless with Black Centre`[iRow],
                           countData$`Colorless with Brown Center`[iRow])

    salmonellaCounts <- rbind(salmonellaCounts, salmonellaCount)
  }

  countData <- cbind(countData, salmonellaCounts)

  return(countData)
}
