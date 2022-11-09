# Function to return the index of the set containing "item" in "overlapSets"
findItem <- function(item, overlapSets) {
	
	# For each set in overlapSets
	for (setIndex in 1:length(overlapSets)) {
		
		# If item is contained within that set, return the index of that set
		if (item %in% overlapSets[[setIndex]]) {
			return(setIndex)
		}
	}
}