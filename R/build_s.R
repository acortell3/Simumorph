
## Function 5. build_s
#' @title Rebuild shape from amplitude and phase values
#' 
#' @description 
#' Reconstructs a shape from given amplitude and phase values of its harmonics. 
#' Optionally returns the Fourier parameters used for reconstruction.
#' 
#' @param x Either a data frame or vector containing amplitude and phase values.  
#'   If `sing.vals = TRUE`, `x` is a data frame with rows equal to the number of harmonics and 4 columns:  
#'   1: Ax (amplitude x), 2: Ay (amplitude y), 3: Phix (phase x), 4: Phiy (phase y).  
#'   If `sing.vals = FALSE`, `x` is a vector of values in the order: Ax1, Ay1, Ax2, Ay2, ..., Axn, Ayn, Phix1, Phiy1, Phix2, Phiy2, ..., Phixn, Phiyn.  
#' @param sing.vals Logical. Indicates if `x` contains singular values as a data frame (`TRUE`) or a combined vector (`FALSE`). Default is `FALSE`.
#' @param fou.pars Logical. If `TRUE`, the function returns a list containing the Fourier parameters and the reconstructed shape coordinates. Default is `FALSE`.
#' @param npts Integer. Number of points to use in the reconstructed shape.
#' 
#' @import Momocs
#' 
#' @return 
#' If `fou.pars = FALSE`, returns a matrix of shape coordinates (x,y).  
#' If `fou.pars = TRUE`, returns a list with two elements:  
#' 1. `Fourier_parameters`: data frame of Fourier coefficients used, with rows for harmonics and columns `an`, `bn`, `cn`, and `dn`.  
#' 2. `Shape`: matrix of shape coordinates (x,y).
#' 
#' @export
build_s <- function(x, sing.vals = F, fou.pars = F, npts){
	
	if (sing.vals == F){ 
		if (isTRUE(is.data.frame(x))){
			if (nrow(x) > 1){
				warning("You've said sing.vals == F, but nrow(x) > 1 (probably you passed an, bn, cn, dn). Are you sure that's what you want?")
			}
			amp_pha <- unlist(x)
		} else {
			amp_pha <- x
		}
		
		## Rebuild parameters from harmonic/phase values
		## Build vectors for selection
		A_vec <- c(1:(length(amp_pha)/2))
		Ph_vec <- c(((length(amp_pha)/2)+1):length(amp_pha))
		
		Ax <- amp_pha[A_vec[A_vec %% 2 == 1]] ## Select odds
		Ay <- amp_pha[A_vec[A_vec %% 2 == 0]] ## Select evens
		Phix <- amp_pha[Ph_vec[Ph_vec %% 2 == 1]]
		Phiy <- amp_pha[Ph_vec[Ph_vec %% 2 == 0]]

		## cos and sin inverted to comply with Momocs
		mat_par_vals <- data.frame("an" = Ax*sin(Phix),
		        		   "bn" = Ax*cos(Phix),
		     		      	   "cn" = Ay*sin(Phiy),
		     	   	      	   "dn" = Ay*cos(Phiy))
	} else {
		if (isTRUE(is.data.frame(x))){
			if (nrow(x) == 1){
				warning("You've said sing.vals == T, but nrow(x) == 1 (probably you passed amplitudes and phases). Are you sure that's what you want?")
			}
		}
		else if (isTRUE(is.vector(x))){
			stop("You've said sing.vals == T, but you provided a vector. Please provide the data in the right format")
		}
		
		Ax <- x[,1]
		Ay <- x[,2]
		Phix <- x[,3]
		Phiy <- x[,4]
		
		## cos and sin inverted to comply with Momocs
		mat_par_vals <- data.frame("an" = Ax*sin(Phix),
		        		   "bn" = Ax*cos(Phix),
		     		      	   "cn" = Ay*sin(Phiy),
		     	   	      	   "dn" = Ay*cos(Phiy))
	}
	
	## Convert to list for momocs
	mat_par_vals_l <- list("an" = mat_par_vals[,1],
			       "bn" = mat_par_vals[,2],
			       "cn" = mat_par_vals[,3],
		               "dn" = mat_par_vals[,4])
	
	# Build shape from pars
	amp_pha_shape <- efourier_i(mat_par_vals_l, nb.h = nrow(mat_par_vals), nb.pts = npts)
	amp_pha_shape <- coo_sample(amp_pha_shape, n = npts)

	if (isTRUE(fou.pars)){
		## Names for pars
		par_names <- rep(NA,nrow(mat_par_vals))
		for (i in 1:nrow(mat_par_vals)){
			par_names[i] <- paste0("harm_",i)
		}
		rownames(mat_par_vals) <- par_names
		res <- list("Fourier_parameters" = mat_par_vals,
			    "Shape" = amp_pha_shape)		
		return(res)
	} else {
		return(amp_pha_shape)
	}

} 

############################################################
## Check that it works
#library(Momocs)

#amp_pha_mat <- readRDS("../../Simu_geos/Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
#shape <- amp_pha_mat[1,]
#shape <- unlist(shape)
#shape <- build_s(shape, sing.vals = F, fou.pars = F, npts = 150) 
#plot(shape[[2]])
#plot(shape)
#shape[[1]]

#str(shape)

#amp_pha_t <- unlist(as.data.frame(amp_pha_mat[1,]))

## Rebuild parameters from harmonic/phase values
## Build vectors for selection
#A_vec <- c(1:(length(amp_pha_t)/2))
#Ph_vec <- c(((length(amp_pha_t)/2)+1):length(amp_pha_t))
		
#Ax <- amp_pha_t[A_vec[A_vec %% 2 == 1]] ## Select odds
#Ay <- amp_pha_t[A_vec[A_vec %% 2 == 0]] ## Select evens
#Phix <- amp_pha_t[Ph_vec[Ph_vec %% 2 == 1]]
#Phiy <- amp_pha_t[Ph_vec[Ph_vec %% 2 == 0]]

#shape_df <- data.frame(Ax,Ay,Phix,Phiy)

#shape_svt <- build_s(shape_df, sing.vals = T, fou.pars = T, npts = 150) 
#plot(shape_svt[[2]])
#shape_svt[[1]]
