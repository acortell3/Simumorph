
## Function 6. proc_dist
#' @title Compute Procrustes distances
#' 
#' @description 
#' Computes Procrustes distances between two shapes after aligning their x and y coordinates using Generalized Procrustes Analysis (GPA). 
#' If `multi = TRUE`, distances are computed between `shape1` and each element in a list of shapes.
#' 
#' @param shape1 A data frame with two columns (x and y coordinates) defining one shape.
#' @param shape2 Either a data frame (if `multi = FALSE`) or a list of data frames (if `multi = TRUE`) with x and y coordinates.
#' @param multi Logical. If `TRUE`, `shape2` must be a list of shapes, and distances are computed from `shape1` to each element. Default is `FALSE`.
#' 
#' @import Momocs
#' @importFrom vegan procrustes
#' 
#' @return 
#' If `multi = FALSE`, returns a numeric value: the Procrustes distance between `shape1` and `shape2`.  
#' If `multi = TRUE`, returns a numeric vector of distances between `shape1` and each shape in the list.  
#' If the list is named, the result will be a named vector.
#' 
#' @export
proc_dist <- function(shape1, shape2, multi=F){
	
	if (multi == F){
		# Center
		shape1_c <- scale(shape1, scale = FALSE, center = TRUE)
		shape2_c <- scale(shape2, scale = FALSE, center = TRUE)
	
		# Scale to unit size
		shape1_sc <- shape1_c / sqrt(sum(shape1_c^2))
		shape2_sc <- shape2_c / sqrt(sum(shape2_c^2))
	
		# Rotate shape2 for best match
		svd_result <- svd(t(shape1_sc) %*% shape2_sc)
		rotation_matrix <- svd_result$v %*% t(svd_result$u)
		shape2_al <- shape2_sc %*% rotation_matrix

		## Compute Procruistes distance
		p_dist <- sqrt(procrustes(shape1_sc,shape2_al,scale = TRUE)$ss)
		names(p_dist) <- "Procrustes distance"
		res <- p_dist
	} else if (multi == T){

		dist_vec <- rep(NA,length(shape2))
		for (i in 1:length(shape2)){
			# Center
			shape1_c <- scale(shape1, scale = FALSE, center = TRUE)
			shape2_c <- scale(as.data.frame(shape2[[i]]), scale = FALSE, center = TRUE)

			# Scale to unit size
			shape1_sc <- shape1_c / sqrt(sum(shape1_c^2))
			shape2_sc <- shape2_c / sqrt(sum(shape2_c^2))
	
			# Rotate shape2 for best match
			svd_result <- svd(t(shape1_sc) %*% shape2_sc)
			rotation_matrix <- svd_result$v %*% t(svd_result$u)
			shape2_al <- shape2_sc %*% rotation_matrix
			
			## Compute Procruistes distance
			dist_vec[i] <- sqrt(procrustes(shape1_sc,shape2_al,scale = TRUE)$ss)
		}
		
		if (!is.null(names(shape2))){
		names(dist_vec) <- paste0("Procrustes distance to ", names(shape2))
		}
		res <- dist_vec
	}
	return(res)
}

## Check that it works
#library(Momocs)
#library(vegan)
#source("build_s.R")
#
#amp_pha_mat <- readRDS("../../Simu_geos/Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
#shape1 <- amp_pha_mat[1,]
#shape2 <- amp_pha_mat[4,]
#shape1 <- build_s(shape1, sing.vals = F, fou.pars = F, npts = 120) 
#shape2 <- build_s(shape2, sing.vals = F, fou.pars = F, npts = 120) 
#
#pdist <- proc_dist(shape1 = shape1, shape2 = shape2, multi = F)
#pdist
#
### Now with multiple distances
#shape1 <- amp_pha_mat[1,]
#shape2 <- amp_pha_mat[2,]
#shape3 <- amp_pha_mat[3,]
#shape4 <- amp_pha_mat[4,]
#shape5 <- amp_pha_mat[5,]
#shape6 <- amp_pha_mat[6,]
#shape7 <- amp_pha_mat[7,]
#shape8 <- amp_pha_mat[8,]
#
#shape1 <- build_s(shape1, sing.vals = F, fou.pars = F, npts = 120) 
#shape2 <- build_s(shape2, sing.vals = F, fou.pars = F, npts = 120) 
#shape3 <- build_s(shape3, sing.vals = F, fou.pars = F, npts = 120) 
#shape4 <- build_s(shape4, sing.vals = F, fou.pars = F, npts = 120) 
#shape5 <- build_s(shape5, sing.vals = F, fou.pars = F, npts = 120) 
#shape6 <- build_s(shape6, sing.vals = F, fou.pars = F, npts = 120) 
#shape7 <- build_s(shape7, sing.vals = F, fou.pars = F, npts = 120) 
#shape8 <- build_s(shape8, sing.vals = F, fou.pars = F, npts = 120) 
#
#shape2 <- list("Shape1" = shape2,
#	       "Shape2" = shape3,
#	       "Shape3" = shape4,
#	       "Shape4" = shape5,
#	       "Shape5" = shape6,
#	       "Shape6" = shape7,
#	       "Shape7" = shape8)
#
#pdist_mult <- proc_dist(shape1 = shape1, shape2 = shape2, multi = TRUE)
#pdist_mult
#
#### And without names
#shape2_nonames <- shape2
#names(shape2_nonames) <- NULL
#pdist_mult_nn <- proc_dist(shape1 = shape1, shape2 = shape2_nonames, multi = TRUE)
#pdist_mult_nn

