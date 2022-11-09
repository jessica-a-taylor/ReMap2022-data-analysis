# Create function that determines whether value a is between values b and c.
betweenFunction <- function(a,b,c) {
	return(b<a & a<c)
}

# Create a function that determines whether two ranges overlap using the between function.
overlapsFunction <- function(S1, E1, S2, E2) {
	if (betweenFunction(S1, S2, E2)) {
		return (TRUE)
	}
	if (betweenFunction(E1, S2, E2)) {
		return (TRUE)
	}
	if (betweenFunction(S2, S1, E1)) {
		return (TRUE)
	}
	if (betweenFunction(E2, S1, E1)) {
		return (TRUE)
	}
	return(FALSE)
}
