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

  for (iRow in 1:nrow(countData))
  {
    eColiCount <- sum(countData$`Red to Pink`[iRow],
                      countData$`Red to Pinkish White`[iRow])

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
    salmonellaCount <- sum(countData$`Partial Green`[iRow], countData$Yellow[iRow])

    if (iRow == 1)
    {
      salmonellaCounts <- salmonellaCount
    }else
    {
      salmonellaCounts <- rbind(salmonellaCounts, salmonellaCount)
    }

  }

  countData <- cbind(countData, salmonellaCounts)

  return(countData)
}
