## Function 3. Create morphospace given a series of shapes
#' @title Create morphospace from shape outlines
#' 
#' @description
#' Creates a morphospace from a set of shapes by extracting their amplitude and phase values.
#' Optionally, it can also compute the covariance matrix of these parameters.
#' 
#' @param x A \code{coo} object (from \pkg{Momocs}) containing the shapes.
#' @param expand Logical. Whether to impute new shapes between existing ones. Default is \code{FALSE}.
#' @param mode Character. Only used if \code{expand = TRUE}. Currently supports \code{"mean"} which imputes mean shapes between pairs of shapes. More modes may be added in the future.
#' @param cov.mat Logical. Whether to compute and return the covariance matrix. Default is \code{TRUE}.
#' @param n.harm Integer. Number of harmonics to use.
#' @param type.names Optional character vector with type names for each shape (length must match \code{length(x)}). If provided, these names will label the rows.
#' 
#' @import Momocs
#' @importFrom matrixcalc is.positive.definite
#' @importFrom stats cov
#' @export
#' 
#' @return
#' If \code{cov.mat = FALSE}, returns a list with:
#' \item{Fourier_parameters}{Fourier parameters of \code{x}.}
#' \item{amp_phase}{Data frame with amplitudes and phases for each shape.}
#' 
#' If \code{cov.mat = TRUE}, returns a list with the above plus:
#' \item{covariance_matrix}{Covariance matrix of the amplitude and phase values.}
morphospace <- function(x,expand = F, mode = "mean", cov.mat = T, n.harm, type.names = NULL){
	## Extract parameters from outlines. To avoid flipping shapes, align before through the x axis (coo_alignxax) and set normalisation to false
	morpho_pars <- efourier(coo_center(coo_scale(coo_alignxax(x))),nb.h = n.harm, norm = F, start = T)

	expanded_ms <- list()
	
	if (expand == T){
		mat_index <- 1 ## Need this to populate the matrix
		for (i in 1:length(x)){
			for (j in i:length(x)){
				shape <- apply(cbind(morpho_pars[i],morpho_pars[j]),1,mean)
				expanded_ms[[mat_index]] <- shape
				mat_index <- mat_index + 1
			}
		}
	} else {
		for (i in 1:length(x)){
			expanded_ms[[i]] <- morpho_pars[i]
		}
	}

	## Extract parameters without constants
	par_vals <- sapply(expanded_ms,function(x) unlist(x))[c(0:n.harm*4),]
	par_vals <- sapply(expanded_ms,function(x) unlist(x))
	
	# Create covariance matrix for amplitude and phase
	# Create vector for column names
	A_names <- c()
	Phi_names <- c()
	for (i in 1:n.harm){
		As <- c(paste0("Ax",i),paste0("Ay",i))
		Phis <- c(paste0("Phix",i),paste0("Phiy",i))

		A_names <- append(A_names,As)
		Phi_names <- append(Phi_names,Phis)
	}

	c_names <- c(A_names,Phi_names)
	
	# Vector for rownames (if types have been specified)
	if (expand == T && !is.null(type.names)){
		## Create rownames for amp_pha_mat
		names_apm <- rep(NA,length(expanded_ms))
		
		index <- 1
		for (i in 1:length(x)){
			for (j in i:length(x)){
				names_apm[index] <- paste0(type.names[i],"_m_",type.names[j])
				index <- index + 1
			}
		}
		type.names <- names_apm
	}

	## Prepare data frame to store amplitude and phase
	amp_pha_mat <- as.data.frame(matrix(0, nrow = length(expanded_ms), ncol = length(c_names), dimnames = list(type.names,c_names)))
n.harm
	## Create list for looping distributing amps and phases
	harms <- vector("list",n.harm)
	for (i in 1:n.harm){
		harms[[i]] <- t(par_vals[c(i,i+n.harm,i+n.harm*2,i+n.harm*3),])
	}	

	## Populate with amplitudes
	index <- 1 ## I need this for correspondence between list and matrix elements
	for (i in 1:length(harms)){
		amp_pha_mat[,index] <- amplitude(harms[[i]], coordinate = "x")
		index <- index + 1
		amp_pha_mat[,index] <- amplitude(harms[[i]], coordinate = "y")
		index <- index + 1	
	}
	
	## Populate with phases
	index <- (ncol(amp_pha_mat)/2)+1 ## I need this for correspondence between list and matrix elements
	for (i in 1:length(harms)){
		amp_pha_mat[,index] <- phase(harms[[i]], coordinate = "x")
		index <- index + 1
		amp_pha_mat[,index] <- phase(harms[[i]], coordinate = "y")
		index <- index + 1	
	}
		
	if (cov.mat == T){
		## Create covariance matrix
		amp_pha_cov <- cov(amp_pha_mat)
		
		if (!isTRUE(is.positive.definite(amp_pha_cov))){
			warning("The covariance matrix is not positive definite. This will likely give you problems further down the line. Consider expanding your morphospace.")
		}

		return(list(morpho_pars,amp_pha_mat,amp_pha_cov))
	} else {
		return(list(morpho_pars,amp_pha_mat))
	}
}


#nh <- 25
#proba <- morphospace(geo_out,expand = T, mode = "mean", cov.mat = T, n.harm = nh, type.names = types)
### Check the morphospace
#amp_pha_mat <- proba[[2]]
#npts <- 120
#n.harm <- nh

#col_vec <- rep("darkslategrey",nrow(amp_pha_mat))
#col_vec[not_selected] <- "khaki3"

## For indexing
#chunks <- list("chunk1" = c(1,100),
#	       "chunk2" = c(101,200),
##	       "chunk3" = c(201,300),
##	       "chunk4" = c(301,400),
##	       "chunk5" = c(401,406))
##index <- 1 # For colouring
#
#for (j in 1:length(chunks)){
#	png(paste0("../Working_figures/Morphospace_",chunks[[j]][1],"_",chunks[[j]][2],".png"), res = 50, height = 1500, width = 1500)
#	par(mfrow = c(10,10))
#	for (i in chunks[[j]][1]:chunks[[j]][2]){
#		## Prepare shape
#		amp_pha <- amp_pha_mat[i,] 
#		amp_pha <- unlist(as.data.frame(amp_pha))
#
#		## Rebuild parameters from harmonic/phase values
#		
#		## Build vectors for selection
#		A_vec <- c(1:(ncol(amp_pha_mat)/2))
#		Ph_vec <- c(((ncol(amp_pha_mat)/2)+1):ncol(amp_pha_mat))
#		
#		Ax <- amp_pha[A_vec[A_vec %% 2 == 1]] ## Select odds
#		Ay <- amp_pha[A_vec[A_vec %% 2 == 0]] ## Select evens
#		Phyx <- amp_pha[Ph_vec[Ph_vec %% 2 == 1]]
#		Phyy <- amp_pha[Ph_vec[Ph_vec %% 2 == 0]]
#
#		## cos and sin inverted to comply with Momocs
#		mat_par_vals <- data.frame("an" = Ax*sin(Phyx),
#		        		   "bn" = Ax*cos(Phyx),
#		     		      	   "cn" = Ay*sin(Phyy),
#		     	   	      	   "dn" = Ay*cos(Phyy))
#		
#		mat_par_vals_l <- list("an" = mat_par_vals[,1],
#			               "bn" = mat_par_vals[,2],
#			               "cn" = mat_par_vals[,3],
#		        	       "dn" = mat_par_vals[,4])
#
#		amp_pha_shape <- efourier_i(mat_par_vals_l, nb.h = n.harm, nb.pts = npts)
#		amp_pha_shape <- coo_sample(amp_pha_shape, n = npts)
#
#		plot(amp_pha_shape, type = "n", asp = 1, xlab = "", ylab = "", main = rownames(amp_pha_mat)[i])
#		polygon(amp_pha_shape, col = col_vec[index], border = "black")	
#		index <- index + 1
#	}
#
#	dev.off()
#}








