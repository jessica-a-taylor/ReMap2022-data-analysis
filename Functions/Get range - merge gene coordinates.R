# Function to merge start and end coordinates

mergeCoordinates <- function(dataToUse) {
  return(paste(dataToUse$start,"-",dataToUse$end, sep = ""))
}
