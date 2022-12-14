# Create function that determines whether value a is between values b and c.
betweenFunction <- function(a,b,c) {
	return(b<a & a<c)
}

# Create a function that determines whether two ranges overlap using the between function.
overlapsFunction <- function(S1, E1, S2, E2) {
  if (betweenFunction(S1, S2, E2)) {
    return (TRUE)
  }
  else if (betweenFunction(E1, S2, E2)) {
    return (TRUE)
  }
  else if (betweenFunction(S2, S1, E1)) {
    return (TRUE)
  }
  else if (betweenFunction(E2, S1, E1)) {
    return (TRUE)
  }
  else return(FALSE)
}

# Create a function that determines the degree to which two ranges overlap using the between function.
newOverlapsFunction <- function(S1, E1, S2, E2) {
  if (betweenFunction(S1, S2, E2)==TRUE & betweenFunction(E1, S2, E2)==TRUE) {
    return(E1-S1)
  }
  else if (betweenFunction(S2, S1, E1) & betweenFunction(E2, S1, E1)) {
    return(E2-S2)
  }
  else if (betweenFunction(S1, S2, E2)) {
    return(E2-S1)
  }
  if (betweenFunction(E1, S2, E2)) {
    return(E1-S2)
  }
  else if (betweenFunction(S2, S1, E1)) {
    return(E1-S2)
  }
  else if (betweenFunction(E2, S1, E1)) {
    return(E2-S1)
  }
  return(0)
}

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

