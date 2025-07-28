
## Function 1. Compute amplitude from an, bn, cn, dn
#' @title Compute amplitude from Fourier parameters
#' 
#' @description 
#' Computes the amplitude from a given set of Fourier parameters \code{an}, \code{bn}, \code{cn}, and \code{dn}.
#' 
#' @param x A data frame or matrix with columns \code{an}, \code{bn}, \code{cn}, and \code{dn}.
#' @param coordinate Character. Specifies whether to compute amplitude for the "x" coordinate, "y" coordinate, or "all". Must be one of \code{"x"}, \code{"y"}, or \code{"all"}.
#' 
#' @return A numeric vector containing amplitude values for each harmonic.
#' 
#' @export
amplitude <- function(x, coordinate=c("x","y","all")){
	coordinate <- match.arg(coordinate)
	if (coordinate == "x"){
		result <- apply(x,1,function(x) sqrt(x[1]^2+x[2]^2))
	} else if (coordinate == "y"){
		result <- apply(x,1,function(x) sqrt(x[3]^2+x[4]^2))
	} else if (coordinate == "all") {
		result <- apply(x,1,function(x) sqrt(x[1]^2+x[2]^2+x[3]^2+x[4]^2))
	}
	return(result)
}
