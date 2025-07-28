
## Function 2. Compute phase from an, bn, cn, dn
#' @title Compute phase from Fourier parameters
#' 
#' @description 
#' Computes the phase from a given set of Fourier parameters \code{an}, \code{bn}, \code{cn}, and \code{dn}.
#' 
#' @param x A data frame or matrix with columns \code{an}, \code{bn}, \code{cn}, and \code{dn}.
#' @param coordinate Character. Specifies whether to compute phase for the \code{"x"} or \code{"y"} coordinate. Must be one of \code{"x"} or \code{"y"}.
#' 
#' @return A numeric vector containing phase values for each harmonic.
#' 
#' @export
phase <- function(x, coordinate = c("x","y")){
	coordinate <- match.arg(coordinate)
	if (coordinate == "x"){
		obj <- x[,c(1,2)]
	} else if (coordinate == "y"){
		obj <- x[,c(3,4)]
	}
	result <- apply(obj,1,function(x) atan2(x[1],x[2]))
	return(result)
}


