
## Function 7. Interpolate shape from two observed shapes
#' @title Interpolate shape given two sets of Fourier parameters
#' 
#' @description 
#' Computes an interpolated shape from a given set of Fourier parameters \code{an}, \code{bn}, \code{cn}, and \code{dn}.
#' 
#' @param shape1 A data frame or matrix with columns \code{an}, \code{bn}, \code{cn}, and \code{dn}. It can also be a list with elements \code{an}, \code{bn}, \code{cn} and \code{dn}
#' @param shape2 A data frame or matrix with columns \code{an}, \code{bn}, \code{cn}, and \code{dn}. It can also be a list with elements \code{an}, \code{bn}, \code{cn} and \code{dn}
#' @param method Character. Specifies whether a single mean shape is computed, a sequence of intermediate shapes between shape1 and shape2, or a weighted mean is computed. More methods will be implemented in the future
#' @param n.shapes. If method = "sequence", the number of interpolated shapes
#' @param fou.pars Logical. If `TRUE`, the function returns a list containing the Fourier parameters and the reconstructed shape coordinates. Default is `FALSE`.
#' @param npts Integer. Number of points to use in the reconstructed shape.
#' @param w If method = "weighted", a vector with two values for the weights for each shape
#'
#' @return 
#' If `fou.pars = FALSE`, returns a matrix of shape coordinates (x,y).  
#' If `fou.pars = TRUE`, returns a list with two elements:  
#' 1. `Fourier_parameters`: data frame of Fourier coefficients used, with rows for harmonics and columns `an`, `bn`, `cn`, and `dn`.  
#' 2. `Shape`: matrix of shape coordinates (x,y).
#' If `method = "sequence"`, returns a list with length = sequence. Each element within the list is structured as above
#' 
#' @export
interpol_s <- function(shape1, shape2, method = c("mean","sequence","weighted"), n.shapes = NULL, npts, fou.pars = F, w = c(0.5,0.5)){
	
	method <- match.arg(method)

	## Check/produce errors/warnings
	if (! (inherits(shape1, "list") || inherits(shape1, "data.frame"))) {
		stop("Error: shape1 must be a list or a data.frame")
	}
	if (! (inherits(shape2, "list") || inherits(shape2, "data.frame"))) {
		stop("Error: shape2 must be a list or a data.frame")
	}
	if (method == "mean" & !is.null(n.shapes)){
		warning("You've chosen method = mean  but n.shapes != NULL. Are you sure that's what you want?")
	} else if (method == "sequence" & is.null(n.shapes)){
		stop("You've chosen method == sequence, but have not provided the number of shapes you want")
	} else if (!(method %in% c("mean","sequence"))){
		stop("Selected method not currently supported")
	}
	
	shape1 <- unlist(shape1)
	shape2 <- unlist(shape2)

	df_ss <- data.frame(shape1,shape2)

	if (method == "mean"){
		## Interpolate value
		new_shape <- apply(df_ss,1,mean)
		n.harm <- length(new_shape)/4 ## Define n.harm from source
		## Assign Fourier parameters
		f.pars <- data.frame("an" = new_shape[1:n.harm],
				     "bn" = new_shape[(1+n.harm):(n.harm*2)],
				     "cn" = new_shape[(1+n.harm*2):(n.harm*3)],
				     "dn" = new_shape[(1+n.harm*3):length(new_shape)])
		## Compute amplitudes and phases and create vector
		amps_x <- amplitude(f.pars, coordinate = "x")
		amps_y <- amplitude(f.pars, coordinate = "y")
		phi_x <- phase(f.pars, coordinate = "x")
		phi_y <- phase(f.pars, coordinate = "y")

		## Vector for amp/phase values
		val_vec <- rep(NA,length(new_shape))
		
		val_vec[seq(1,length(val_vec)-1, by = 2)] <- c(amps_x,phi_x)
		val_vec[seq(2,length(val_vec), by = 2)] <- c(amps_y,phi_y)
		
		new_shape <- build_s(val_vec, sing.vals = F, fou.pars = fou.pars, npts = npts)
		return(new_shape)
	} else if (method == "sequence"){
		new_shape_df <- apply(df_ss,1,function(x) seq(x[1], x[2],length.out = n.shapes))
		n.harm <- ncol(new_shape_df)/4 ## Define n.harm from source

		## create list to store values
		new_shapes <- list()
		
		## Populate list with interpolated shapes
		for (i in 1:n.shapes){
			## Interpolate value
		new_shape <- new_shape_df[i,]
			## Assign Fourier parameters
			f.pars <- data.frame("an" = new_shape[1:n.harm],
					     "bn" = new_shape[(1+n.harm):(n.harm*2)],
				       	     "cn" = new_shape[(1+n.harm*2):(n.harm*3)],
				             "dn" = new_shape[(1+n.harm*3):length(new_shape)])
			
			## Compute amplitudes and phases and create vector
			amps_x <- amplitude(f.pars, coordinate = "x")
			amps_y <- amplitude(f.pars, coordinate = "y")
			phi_x <- phase(f.pars, coordinate = "x")
			phi_y <- phase(f.pars, coordinate = "y")

			## Vector  amp/phase values
			val_vec <- rep(NA,length(new_shape))
		
			val_vec[seq(1,length(val_vec)-1, by = 2)] <- c(amps_x,phi_x)
			val_vec[seq(2,length(val_vec), by = 2)] <- c(amps_y,phi_y)

			new_shapes[[i]] <- build_s(val_vec, sing.vals = F, fou.pars = fou.pars, npts = npts)
		} 
		return(new_shapes)
	} else if (method == "weighted"){
		
		## Interpolate value
		new_shape <- apply(df_ss,1,function(x) sum(x*w) / sum(w))
		n.harm <- length(new_shape)/4 ## Define n.harm from source
		## Assign Fourier parameters
		f.pars <- data.frame("an" = new_shape[1:n.harm],
				     "bn" = new_shape[(1+n.harm):(n.harm*2)],
				     "cn" = new_shape[(1+n.harm*2):(n.harm*3)],
				     "dn" = new_shape[(1+n.harm*3):length(new_shape)])
		## Compute amplitudes and phases and create vector
		amps_x <- amplitude(f.pars, coordinate = "x")
		amps_y <- amplitude(f.pars, coordinate = "y")
		phi_x <- phase(f.pars, coordinate = "x")
		phi_y <- phase(f.pars, coordinate = "y")

		## Vector for amp/phase values
		val_vec <- rep(NA,length(new_shape))
		
		val_vec[seq(1,length(val_vec)-1, by = 2)] <- c(amps_x,phi_x)
		val_vec[seq(2,length(val_vec), by = 2)] <- c(amps_y,phi_y)
		
		new_shape <- build_s(val_vec, sing.vals = F, fou.pars = fou.pars, npts = npts)
		return(new_shape)
	}
	else if (!(method %in% c("mean","sequence","weighted"))){
		stop("Selected method not currently supported")
	}
}

### Check that it works
#library(Momocs)
#library(Simumorph)
#mpars <- readRDS("../../Simu_geos/Utilities/geo_out.rds")

#geo_out_li <- slice(mpars, c(26,27))
#morpho_pars_li <- list()
#morpho_pars_li <- efourier(coo_center(coo_scale(coo_alignxax(geo_out_li))),nb.h = n.harm, norm = F, start = T)

## Interpolate linear values for three middle shapes

#s1 <- data.frame("an" =  morpho_pars_li[1][1:7],
#		 "bn" =  morpho_pars_li[1][8:14],
#		 "cn" =  morpho_pars_li[1][15:21],
#		 "dn" =  morpho_pars_li[1][22:28])
#s2 <- list("an" =  morpho_pars_li[2][1:7],
#	   "bn" =  morpho_pars_li[2][8:14],
#	   "cn" =  morpho_pars_li[2][15:21],
#	   "dn" =  morpho_pars_li[2][22:28])

## Try with signgle shape
#int_shape <- interpol_s(shape1 = s1, shape2 = s2, method = "maen", npts = 120, fou.pars = T)
#plot(int_shape[[2]], type = "l")

## Try with sequence
#int_shape <- interpol_s(shape1 = s1, shape2 = s2, method = "sequence", npts = 120, fou.pars = T, n.shapes = 3)

#par(mfrow = c(1,3))
#plot(int_shape[[1]]$Shape, type = "l", asp = 1)
#plot(int_shape[[2]]$Shape, type = "l", asp = 1)
#plot(int_shape[[3]]$Shape, type = "l", asp = 1)





