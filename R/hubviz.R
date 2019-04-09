#'
#' hubViz: A novel tool for Hub-centric Visualization
#' 
#' This package provides funtions for fitting hubviz, a novel tool for hub-centric visualization
#' 
#' @export 
#' 
#' @param dataset Matrix of binary values with n samples and p variables
#' @param nsample Number of samples to analyze (default: use all samples)
#' @param nitem Number of variables to analyze (default: use all variables)
#' @param ndim Dimension of latent position (default: 2)
#' @param niter Total number of iterations in MCMC (default: 30000)
#' @param nburn Number of burn-in iterations in MCMC (default: 5000)
#' @param nthin Thinning number of MCMC (default: 5)
#' @param nprint Trace size of MCMC (default: 100)
#' @param jump_theta Jump size for the intercept parameter (default: 1.0)
#' @param jump_w Jump size for latent space (default: 0.2)
#' @param pr_mean_theta Prior mean for the intercept parameter (default: 0.0)
#' @param pr_sd_theta Prior Standard deviation for the intercept parameter (default: 10.0)
#' @param pr_mean_w Prior mean for latent space (default: 0.0)
#' @param prior_a Prior Shape parameter for variance of latent space (default: 0.001)
#' @param prior_b Prior Scale parameter for variance of latent space (default: 0.001)
#' @param option Logical. If TRUE then print MCMC iteration (default: FALSE)
#' @param verbose Logical. If TRUE then return MCMC results (default: FALSE)
#' @param cores Number of CPUs to be used for computing (default: 1)
#' 
#' @examples
#' n <- 100
#' p <- 10
#' bmat <- matrix( 0, n, p )
#' for ( j in 1:ncol(bmat) ) {
#'  bmat[ 1:(10*j), j ] <- 1
#' }
#' rownames( bmat ) <- paste( "Sample", 1:n, sep="" )
#' colnames( bmat ) <- paste( "V", 1:p, sep="" )
#' 
#' # Fit hubviz
#' nsamp <- nrow( bmat ); nitem <- ncol( bmat )
#' fit <- hubviz( bmat, nsample = nsamp, nitem = nitem)
#' fit
#' 
#' # Extraction of results
#' estimate( fit )
#' 
#' # Hub-centric plot
#' plot( fit )


hubviz <- function(dataset, nsample=nrow(dataset), nitem=ncol(dataset), ndim = 2, niter = 30000, nburn = 5000, nthin = 5, nprint = 100,
                    jump_theta = 1.0, jump_w = 0.2, pr_mean_theta = 0.0, pr_sd_theta = 10.0, pr_mean_w = 0.0, 
                    prior_a = 0.001, prior_b = 0.001, option = FALSE, verbose=FALSE, cores = 1){

  message("\n")
  message("---------------------\n")
  message("Running MCMC... \n")
  message("---------------------\n")
  
  if((niter - nburn) %% nthin == 0){

        output = hubviz_cpp(dataset, nsample, nitem, ndim, niter, nburn, nthin, nprint,
                      jump_theta, jump_w, pr_mean_theta, pr_sd_theta, pr_mean_w, 
                      prior_a, prior_b, option, cores)
    
    nmcmc = as.integer((niter - nburn) / nthin)
    max.address = which.max(output$posterior)
    w.star = output$w[max.address,,]
    w.proc = array(0,dim=c(nmcmc,nitem,ndim))
    
    for(iter in 1:nmcmc){
      w.iter = output$w[iter,,]
      if(iter != max.address) w.proc[iter,,] = procrustes(w.iter,w.star)$X.new
      else w.proc[iter,,] = w.iter
    }
    
    theta.est = apply(output$theta,2,mean)
    
    w.est = matrix(NA,nitem,ndim)
    for(i in 1:nitem){
      for(j in 1:ndim){
        w.est[i,j] = mean(w.proc[,i,j])
      }
    }
    sigma.w = mean(output$sigma_w)
    
    result = list(theta=output$theta, w=w.proc, 
                           theta.estimate=theta.est, w.estimate=w.est, 
                           sigma.w = output$sigma_w, sigma.w.estimate = sigma.w,
                           accept_theta=output$accept_theta,
                           accept_w=output$accept_w)

  }
  else{
    print("Error: The total size of MCMC sample is not integer")
    return(-999)
  }

   if (verbose) {
   result <- result
   } else {
   result <- list( theta.estimate=theta.est, w.estimate=w.est, sigma.w.estimate = sigma.w )
   }

	methods::new("hubviz",
	             data = dataset,
	             init = list(nsample=nsample, nitem=nitem, ndim = ndim, niter = niter, nburn = nburn, nthin=nthin),
	             result = result
	             )
  }
  