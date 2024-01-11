# numericCounts
#
# Convert plate counts to numeric form
#
# file numericCounts
#
# inputs
# 	countData     - Initial Data to be converted

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

numericCounts <- function(countData) {

  # Ensure all the counts data is numeric
  if (!is_empty(which(colnames(countData) == "Black/Brown")))
  {
    countData$`Black/Brown` <- as.numeric(countData$`Black/Brown`)
  }

  countData$Colorless <- as.numeric(countData$Colorless)

  if (!is_empty(which(colnames(countData) == "Dark Blue")))
  {
    countData$`Dark Blue` <- as.numeric(countData$`Dark Blue`)
  }

  if (!is_empty(which(colnames(countData) == "Dark Pink")))
  {
    countData$`Dark Pink` <- as.numeric(countData$`Dark Pink`)
  }

  if (!is_empty(which(colnames(countData) == "Light Pink")))
  {
    countData$`Light Pink` <- as.numeric(countData$`Light Pink`)
  }

  if (!is_empty(which(colnames(countData) == "Partial Green")))
  {
    countData$`Partial Green` <- as.numeric(countData$`Partial Green`)
  }

  countData$`Pink, Bile precipitate` <- as.numeric(countData$`Pink, Bile precipitate`)

  if (!is_empty(which(colnames(countData) == "Purple/Voilet")))
  {
    countData$`Purple/Voilet` <- as.numeric(countData$`Purple/Voilet`)
  }

  if (!is_empty(which(colnames(countData) == "Salmon to Red")))
  {
    countData$`Salmon to Red` <- as.numeric(countData$`Salmon to Red`)
  }

  if (!is_empty(which(colnames(countData) == "Turquoise")))
  {
    countData$`Turquoise` <- as.numeric(countData$`Turquoise`)
  }

  if (!is_empty(which(colnames(countData) == "Turquoise/Light Blue")))
  {
    countData$`Turquoise/Light Blue` <- as.numeric(countData$`Turquoise/Light Blue`)
  }

  if (!is_empty(which(colnames(countData) == "Violet/Purple")))
  {
    countData$`Violet/Purple` <- as.numeric(countData$`Violet/Purple`)
  }

  countData$Yellow <- as.numeric(countData$Yellow)
  countData$`Colorless with Black center` <- as.numeric(countData$`Colorless with Black center`)
  countData$`Colorless with Brown Center` <- as.numeric(countData$`Colorless with Brown Center`)
  countData$`Colorless with Yellow center` <- as.numeric(countData$`Colorless with Yellow center`)
  countData$White <- as.numeric(countData$White)

  countData$Total <- as.numeric(countData$Total)

  return(countData)
}
