
## Function 8. Check the number of necessary shapes to produce a positive definite matrix
#' @title Check matrix validity
#' 
#' @description
#' Given a full covariance matrix, it assesses how many samples are necessary, so that it is positive definite
#' 
#' @param x The matrix with the full morphospace
#' @param k the number of eigen values to include. Must be higher than 1
#' @param rnd Number of randomisations. Default is 100
#' @param plotvalues Boolean. Whether the user wants to visualise also the eigenvalues. If F, it only plots the dot product of the eigenvectors and the minimum eigenvalues. Default is F
#'
#' @export
#' 
#' @return
#' A plot showing the directions of the eigenvectors and the minimum eigenvalues. If \code{plotvalues = TRUE} it also shows the k eigenvalues 

valid_morphospace <- function(x,k,rnd,plotvalues = T){

	if (k == 1){
		stop("k must be higher than 1")
	}
	
	## Actual function
	vars <- ncol(cov(x))
	target <- mat[sample(1:nrow(x),vars),]
	cov_target <- cov(target)
	eigen_target <- eigen(cov_target)$vectors[,1:k]
	
	eigen_vals <- data.frame(replicate(k+3,numeric()))
	colnames(eigen_vals) <- c("eVector",paste0("eValue",c(1:k)),"eValuemin","Obs")
	
	for (h in 1:rnd){
		eigen_vals_prov <- data.frame(replicate(k+3,numeric()))
		colnames(eigen_vals_prov) <- c("eVector",paste0("eValue",c(1:k)),"eValuemin","Obs")
		
		for (i in vars:nrow(x)){
			sp_mat <- x[sample(1:nrow(x),i),]
			cvmat <- cov(sp_mat)
			eigen_cvmat <- eigen(cvmat)
			
			# Compute eigen vectors	
			eigen_vals_prov[i-27,1] <- mean(abs(colSums(eigen_target[,1:k]*eigen_cvmat$vector[,1:k])))
		
			# Distribute eigen values
			eigen_vals_prov[i-27,2:(k+1)] <- eigen_cvmat$values[1:k]
			eigen_vals_prov[i-27,k+2] <- min(eigen_cvmat$values)
		
			#Observation number	
			eigen_vals_prov[i-27,k+3] <- i
		}
		eigen_vals <- rbind(eigen_vals,eigen_vals_prov)
	}
	
	if (isTRUE(plotvalues)){
		
		if (k > 7){
			warning("k is too high to fit into a single plot. You don't need that many eigenvalues anyway...")
		}
		if (k <= 4){
			par(mfrow = c(round(k/2)+1,2))
		} else {
			par(mfrow = c(round(sqrt(k+2)),round(sqrt(k+2))))
		}
		plot(x = eigen_vals$Obs, y = eigen_vals$eVector, pch = 16, col = adjustcolor("hotpink3", alpha.f = 0.2), xlab = "n samples", ylab = "value", main = paste0("Eigenvector dot product (k = ",k,")"))
		plot(x = eigen_vals$Obs, y = eigen_vals$eValuemin, col = adjustcolor("olivedrab4", alpha.f = 0.2), xlab = "n samples", ylab = "value", pch = 16, main = "Minimum eigenvalue")
		polygon(x = c(0,0,nrow(mat)+50,nrow(mat)+50), y = c(-0.05,1e-6,1e-6,-0.05), col = adjustcolor("pink",alpha.f=0.2), border = "red", lty =2)
		for (i in 1:k){
			plot(x = eigen_vals$Obs, y = eigen_vals[,i+1], pch = 16, col = adjustcolor("lightseagreen", alpha.f = 0.2), xlab = "n samples", ylab = "value", main = paste0("eigenvalue",i))
		}
	} else {
		par(mfrow = c(1,2))
		plot(x = eigen_vals$Obs, y = eigen_vals$eVector, pch = 16, col = adjustcolor("hotpink3", alpha.f = 0.2), xlab = "n samples", ylab = "value", main = paste0("Eigenvector dot product (k = ",k,")"))
		plot(x = eigen_vals$Obs, y = eigen_vals$eValuemin, col = adjustcolor("olivedrab4", alpha.f = 0.2), xlab = "n samples", ylab = "value", pch = 16, main = "Minimum eigenvalue")
		polygon(x = c(0,0,nrow(mat)+50,nrow(mat)+50), y = c(-0.05,1e-6,1e-6,-0.05), col = adjustcolor("pink",alpha.f=0.2), border = "red", lty =2)
	}
}
	
#valid_morphospace(x = mat, k =1 , rnd = rnd,plotvalues = T)
	
	






