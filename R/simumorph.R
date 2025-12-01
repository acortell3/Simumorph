

## Function 4. simumorph
#' @title Simulate shape evolution in morphospace given a covariance matrix
#' 
#' @description
#' Simulates shapes over \code{sim} iterations based on a covariance matrix derived from a morphospace.
#' The simulation starts from an initial shape and optionally converges or varies relative to target shapes.
#' 
#' The morphospace includes amplitude and phase values of shapes. Users can add shapes not originally in the morphospace by providing their amplitude and phase values.
#' Note that added shapes are excluded from the covariance matrix unless explicitly included before covariance calculation.
#' 
#' @param x A positive definite covariance matrix of shape parameters from the morphospace.
#' @param m.space A data frame of amplitude and phase values for all shapes in the morphospace (ideally from \code{morphospace()}).
#' @param init Integer or character. The initial shape's row index or name in \code{m.space} where simulation begins.
#' @param init.from.morphospace Logical. Whether to include the initial shape in truncations. Default is \code{TRUE}.
#' @param target Integer, character, or vector thereof. Row index(es) or name(s) in \code{m.space} of target shape(s) guiding the simulation.
#' @param target.from.morphospace Logical. Whether to include target shape(s) in truncations. Default is \code{TRUE}.
#' @param method Character. Simulation mode; one of \code{"AtoA"}, \code{"AtoB"}, \code{"AtoMult"}, or \code{"Free"}. Defaults to \code{"AtoA"}.
#' \describe{
#'   \item{\code{"AtoA"}}{Simulation stays near the initial shape.}
#'   \item{\code{"AtoB"}}{Simulation converges towards the target shape(s).}
#'   \item{\code{"AtoMult"}}{Simulation allows variation but remains close to multiple targets (e.g., full morphospace).}
#'   \item{\code{"Free"}}{Simulation is unconstrained and may deviate freely from known shapes.}
#' }
#' @param sim Integer. Number of simulation iterations.
#' @param npts Integer. Number of points to reconstruct each shape.
#' @param a Numeric. Phase local threshold to smooth changes. Default is 0.2.
#' @param e Numeric. Procrustes distance threshold for acceptance. Default is 0.15.
#' @param dynamic_e Logical or data frame. If \code{TRUE}, modifies distance threshold \code{e} dynamically during simulation based on a two-column data frame with time steps and threshold values. Default is \code{FALSE}.
#' @param f Numeric. Factor dividing the covariance matrix.
#' @param int.allowed Logical. Whether to allow line intersections in candidate shapes. Default is \code{FALSE}.
#' @param only.shapes Logical. If \code{TRUE}, returns only simulated shapes at each time step. If \code{FALSE}, returns detailed info including parameters and distances. Default is \code{FALSE}.
#' @param max.attempts Integer. Maximum attempts allowed without success in meeting \code{e} threshold before simulation stops. Default is 500.
#' 
#' @import sf
#' @importFrom tmvtnorm rtmvnorm
#' @export
#' 
#' @return
#' A list containing simulated shapes (and optionally parameters and distance metrics) for each iteration.
simumorph <- function(x, m.space, init, init.from.morphospace = T, target, target.from.morphospace = T, method = c("AtoA","AtoB","AtoMult","Free"), sim, npts, a = 0.2, e = 0.15, dynamic_e = F, f = 100, int.allowed = F, only.shapes = F, max.attempts = 500){

	## Define starting shape
	i_shape <- m.space[init,]

	## Define target(s)
	## For the case of the free simulation
	if (method == "Free" | method == "AtoA"){
		target <- init		
		tar_vals <- m.space[init,]
		tar_shape <- build_s(unlist(m.space[init,]), fou.pars = F, npts = npts)
	} else if (method == "AtoB"){
		tar_vals <- m.space[target,]
		tar_shape <- build_s(unlist(m.space[target,]), fou.pars = F, npts = npts)
	} else if (method == "AtoMult"){
		tar_shape <- list()
		for (i in 1:length(target)){
			tar_shape[[i]] <- build_s(unlist(m.space[target[i],]), fou.pars = F, npts = npts)
		}
	} else {
		stop("No valid method provided")
	} 

	## Remove init from morphospace, if so specified
	if (init.from.morphospace == F){
		m.space <- m.space[-init,]
	}

	## Remove target(s) from morphospace, if so specified
	if (target.from.morphospace == F){
		m.space <- m.space[-target,]
	}
	
	## Define truncations for multivariate normals
	# Thresholds for truncated multivariate normals
	upper_thres <- apply(m.space,2,max)
	lower_thres <- apply(m.space,2,min)
	
	## Retrieve number of harmonics
	n.harm <- ncol(m.space)/4
		
	## Build objects to store values
	## Amplitudes and phases
	Axmat <- matrix(nrow = n.harm, ncol = sim)
	rownames(Axmat) <- paste0("harm_",seq(1:n.harm))
	colnames(Axmat) <- paste0("Ax_t_",seq(1:sim)) 
	Aymat <- matrix(nrow = n.harm, ncol = sim)
	rownames(Aymat) <- paste0("harm_",seq(1:n.harm))
	colnames(Aymat) <- paste0("Ay_t_",seq(1:sim)) 
	Phixmat <- matrix(nrow = n.harm, ncol = sim)
	rownames(Phixmat) <- paste0("harm_",seq(1:n.harm))
	colnames(Phixmat) <- paste0("Phix_t_",seq(1:sim)) 
	Phiymat <- matrix(nrow = n.harm, ncol = sim)
	rownames(Phiymat) <- paste0("harm_",seq(1:n.harm))
	colnames(Phiymat) <- paste0("Phiy_t_",seq(1:sim)) 

	## Distances
	proc_distances <- rep(NA,sim)
	names(proc_distances) <- paste0("P.dist_t_",seq(1:sim))
	proc_distances_mult <- matrix(nrow = length(target), ncol = sim) ## In case A to Mult
	rownames(proc_distances_mult) <- paste0("target_",seq(1:length(target)))
	colnames(proc_distances_mult) <- paste0("P.dist_t_",seq(1:sim))

	## Fourier parameters
	par_sim <- list()

	## Store shapes
	simuls <- list()

	## Start actual simulation
	index <- 1
	attempt <- 0 ## To stop the simulation in case several unsuccessful candidates have been proposed
	while (index <= sim){
		## Set e if it is dynamic
		if (!isFALSE(dynamic_e)){
			match_row <- which(dynamic_e[,1] == index)
			if (length(match_row) > 0){
				e <- dynamic_e[match_row,2]
			}
		}

		## Set lower and upper threshold of the phases with a small kernel from the observed phases to avoid large jumps
		lower_thres[((ncol(m.space)/2)+1):ncol(m.space)] <- unlist(as.data.frame(i_shape[((ncol(m.space)/2)+1):ncol(m.space)])) - a
		upper_thres[((ncol(m.space)/2)+1):ncol(m.space)] <- unlist(as.data.frame(i_shape[((ncol(m.space)/2)+1):ncol(m.space)])) + a
 	
		## Bound to [-pi,pi]
		lower_thres[((ncol(m.space)/2)+1):ncol(m.space)] <- pmin(pmax(lower_thres[((ncol(m.space)/2)+1):ncol(m.space)],-pi),pi)
		upper_thres[((ncol(m.space)/2)+1):ncol(m.space)] <- pmin(pmax(upper_thres[((ncol(m.space)/2)+1):ncol(m.space)],-pi),pi)
		
		## Propose candidate values
		cand_vals <- rtmvnorm(1,mean = unlist(as.data.frame(i_shape)), sigma = x/f, lower = lower_thres, upper = upper_thres)

		## Build candidate shape
		cand_shape <- build_s(cand_vals, fou.pars = T, npts = npts)
		closed_shape <- rbind(cand_shape[[2]],cand_shape[[2]][1,]) ## Close for spatial validity and comparison
		closed_shape <- st_polygon(list(closed_shape))
		cross <- st_is_valid(closed_shape) ## Are there interesections

		if (int.allowed == F){
			if (isTRUE(cross)){
				if (method == "AtoA"){
					## Compute procrustes distance 
					pdist <- proc_dist(cand_shape[[2]],tar_shape,multi = F) 
					
					if (pdist <= e){
						## Store values for posterior plotting and assessment
						# Vectors to store amplitude and phase values
						amplitudes <- cand_vals[1:(length(cand_vals)/2)]
						phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
						simuls[[index]] <- cand_shape[[2]]
						Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
						Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
						Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
						Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
						proc_distances[index] <- pdist
						par_sim[[index]] <- cand_shape[[1]]
						i_shape <- cand_vals
			
						print(index) ## To see progress
						index <- index + 1
						attempt <- 0
					} else {
						attempt <- attempt + 1
						if (attempt >= max.attempts){
							stop(paste0("Simulation failed after ", max.attempts," attempts. Consider increasing e"))
						}
					}
				} else if (method == "AtoB"){
					old_e <- e
					## e procrustes distance 
					pdist <- proc_dist(cand_shape[[2]],tar_shape,multi = F) 
					
					if (pdist <= old_e){
						## Store values for posterior plotting and assessment
						# Vectors to store amplitude and phase values
						amplitudes <- cand_vals[1:(length(cand_vals)/2)]
						phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
						simuls[[index]] <- cand_shape[[2]]
						Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
						Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
						Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
						Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
						proc_distances[index] <- pdist
						par_sim[[index]] <- cand_shape[[1]]
						i_shape <- cand_vals
			
						e <- pdist ## Update distance
						print(index) ## To see progress
						index <- index + 1
						attempt <- 0
					} else {
						attempt <- attempt + 1
						if (attempt >= max.attempts){
							e <- e+e*0.1
							warning(paste0("After ", max.attempts," attempts at t = ",index," the simulation was stuck, and e was increased by e*0.1"))
						}
					}
				} else if (method == "AtoMult"){
					## Compute procrustes distance 
					pdist <- proc_dist(cand_shape[[2]], tar_shape, multi = T)

					if (any(pdist <= e)){
						## Store values for posterior plotting and assessment
						# Vectors to store amplitude and phase values
						amplitudes <- cand_vals[1:(length(cand_vals)/2)]
						phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
						simuls[[index]] <- cand_shape[[2]]
						Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
						Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
						Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
						Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
						proc_distances_mult[index] <- pdist
						par_sim[[index]] <- cand_shape[[1]]
						i_shape <- cand_vals
			
						print(index) ## To see progress
						index <- index + 1
						attempt <- 0
					} else {
						attempt <- attempt + 1
						if (attempt >= max.attempts){
							stop(paste0("Simulation failed after ", max.attempts," attempts. Consider increasing e"))
						}
					}
				} else if (method == "Free"){
					## Compute procrustes distance 
					pdist <- proc_dist(cand_shape[[2]],tar_shape,multi = F) 
					
					## Store values for posterior plotting and assessment
					# Vectors to store amplitude and phase values
					amplitudes <- cand_vals[1:(length(cand_vals)/2)]
					phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
					simuls[[index]] <- cand_shape[[2]]
					Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
					Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
					Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
					Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
					proc_distances[index] <- pdist
					par_sim[[index]] <- cand_shape[[1]]
					i_shape <- cand_vals
		
					print(index) ## To see progress
					index <- index + 1
		
				}
			}
		} else if (int.allowed == T){
			if (method == "AtoA"){
				## Compute procrustes distance 
				pdist <- proc_dist(cand_shape[[2]],tar_shape,multi = F) 
					
				if (pdist <= e){
					## Store values for posterior plotting and assessment
					# Vectors to store amplitude and phase values
					amplitudes <- cand_vals[1:(length(cand_vals)/2)]
					phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
					simuls[[index]] <- cand_shape[[2]]
					Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
					Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
					Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
					Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
					proc_distances[index] <- pdist
					par_sim[[index]] <- cand_shape[[1]]
					i_shape <- cand_vals
			
					print(index) ## To see progress
					index <- index + 1
					attempt <- 0
				} else {
					attempt <- attempt + 1
					if (attempt >= max.attempts){
						stop(paste0("Simulation failed after ", max.attempts," attempts. Consider increasing e"))
					}
				}
			} else if (method == "AtoB"){
				old_e <- e
				## e procrustes distance 
				pdist <- proc_dist(cand_shape[[2]],tar_shape,multi = F) 
				
				if (pdist <= old_e){
					## Store values for posterior plotting and assessment
					# Vectors to store amplitude and phase values
					amplitudes <- cand_vals[1:(length(cand_vals)/2)]
					phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
					simuls[[index]] <- cand_shape[[2]]
					Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
					Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
					Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
					Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
					proc_distances[index] <- pdist
					par_sim[[index]] <- cand_shape[[1]]
					i_shape <- cand_vals
		
					e <- pdist ## Update distance
					print(index) ## To see progress
					index <- index + 1
					attempt <- 0
				} else {
					attempt <- attempt + 1
					if (attempt >= max.attempts){
						e <- e+e*0.1
						warning(paste0("After ", max.attempts," attempts at t = ",index," the simulation was stuck, and e was increased by e*0.1"))
					}
				}
			} else if (method == "AtoMult"){
				## Compute procrustes distance 
				pdist <- proc_dist(cand_shape[[2]], tar_shape, multi = T)
				if (any(pdist <= e)){
					## Store values for posterior plotting and assessment
					# Vectors to store amplitude and phase values
					amplitudes <- cand_vals[1:(length(cand_vals)/2)]
					phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
					simuls[[index]] <- cand_shape[[2]]
					Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
					Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
					Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
					Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
					proc_distances_mult[,index] <- pdist
					par_sim[[index]] <- cand_shape[[1]]
					i_shape <- cand_vals
		
					print(index) ## To see progress
					index <- index + 1
					attempt <- 0
				} else {
					attempt <- attempt + 1
					if (attempt >= max.attempts){
						stop(paste0("Simulation failed after ", max.attempts," attempts. Consider increasing e"))
					}
				}
			} else if (method == "Free"){
				## Compute procrustes distance 
				pdist <- proc_dist(cand_shape[[2]],tar_shape,multi = F) 
				
				## Store values for posterior plotting and assessment
				# Vectors to store amplitude and phase values
				amplitudes <- cand_vals[1:(length(cand_vals)/2)]
				phases <- cand_vals[((length(cand_vals)/2)+1):length(cand_vals)]
				simuls[[index]] <- cand_shape[[2]]
				Axmat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 1] 
				Aymat[,index] <- amplitudes[(1:length(amplitudes)) %% 2 == 0]
				Phixmat[,index] <- phases[(1:length(phases)) %% 2 == 1]
				Phiymat[,index] <- phases[(1:length(phases)) %% 2 == 0]
				proc_distances[index] <- pdist
				par_sim[[index]] <- cand_shape[[1]]
				i_shape <- cand_vals
	
				print(index) ## To see progress
				index <- index + 1
	
			}
		}
	}
	if (isTRUE(only.shapes)){
		names(simuls) <- paste0("Shape_t_",seq(1:length(simuls)))
		return(simuls)
	} else {
		names(simuls) <- paste0("Shape_t_",seq(1:length(simuls)))
		names(par_sim) <- paste0("Shape_t_",seq(1:length(par_sim)))
		res <- list("Shapes" = simuls,
			    "Xamplitudes" = Axmat,
			    "Yamplitudes" = Aymat,
			    "Xphases" = Phixmat,
			    "Yphases" = Phiymat,
			    "P.distances" = proc_distances,
			    "F.parameters" = par_sim)

		return(res)
	}
}

####### Check that it works

######################################################
#### Load libraries and utilities
#######################################################

## Libraries
# library(Momocs) ## For GMM functions
# library(vegan) ## For manual Procrustes
# library(tmvtnorm) ## For truncated multivariate normal
# library(sf) ## to check internal intersections
# 
# ## Load necessary functions
# lapply(c("amplitude.R","phase.R","proc_dist.R","build_s.R","morphospace.R"), source)
# 
# 
# ## Load necessary objects (produced with morphospace.R)
# geo_out <- readRDS("../../Simu_geos/Utilities/geo_out.rds") ## Observed shapes
# morpho_pars <- readRDS("../../Simu_geos/Utilities/morpho_pars.rds") ## Observed parameters
# amp_pha_mat <- readRDS("../../Simu_geos/Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
# amp_pha_cov <- readRDS("../../Simu_geos/Utilities/amp_pha_cov.rds") ## Covariance matrix
# 
# 
# ## Examples
# set.seed(1)
# sims <- 100
# 
# out <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "AtoA", sim = sims, npts = 120, only.shapes = T, a = 0.2, e = 0.05, max.attempts = 500, f = 100)
# 
# mat <- matrix(c(seq(1,100),rep(101,20)), nrow = 12, byrow = T) 
# 
# png("../Working_images/AtoA.png", res = 50, height = 1500, width = 1500)
# layout(mat)
# for (i in 1:sims){
# 	plot(out[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
# 	polygon(out[[i]], col = "seagreen")
# }
# dev.off()
# 
# out <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "AtoB", sim = sims, npts = 120, only.shapes = T, a = 0.5, e = 0.5, max.attempts = 500, f = 100)
# 
# mat <- matrix(c(seq(1,100),rep(101,20)), nrow = 12, byrow = T) 
# 
# png("../Working_images/AtoB.png", res = 50, height = 1500, width = 1500)
# layout(mat)
# for (i in 1:sims){
# 	plot(out[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
# 	polygon(out[[i]], col = "seagreen")
# }
# dev.off()
# 
# out <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = c(1,10,50,80), method = "AtoMult", sim = sims, npts = 120, only.shapes = T, a = 1, e = 0.1, max.attempts = 500, f = 100)
# 
# mat <- matrix(c(seq(1,100),rep(101,20)), nrow = 12, byrow = T) 
# 
# png("../Working_images/AtoMult.png", res = 50, height = 1500, width = 1500)
# layout(mat)
# for (i in 1:sims){
# 	plot(out[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
# 	polygon(out[[i]], col = "seagreen")
# }
# dev.off()
# 
# out <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "Free", sim = sims, npts = 120, only.shapes = T, a = 0.2, e = 0.05, max.attempts = 500, f = 100)
# 
# mat <- matrix(c(seq(1,100),rep(101,20)), nrow = 12, byrow = T) 
# 
# png("../Working_images/Free.png", res = 50, height = 1500, width = 1500)
# layout(mat)
# for (i in 1:sims){
# 	plot(out[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
# 	polygon(out[[i]], col = "seagreen")
# }
# dev.off()
# 
# ## Try dynamic e
# out_de <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, method = "AtoA", sim = sims, npts = 120, only.shapes = T, e = 1, dynamic_e = data.frame("time" = c(3,6,10),
# 																			      "e" = c(1,0.5,0.4)))
# 
# png("../Working_images/dynamic_e.png", res = 50, height = 1500, width = 1500)
# layout(mat)
# for (i in 1:sims){
# 	plot(out_de[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
# 	polygon(out_de[[i]], col = "seagreen")
# }
# dev.off()
# 
# 
