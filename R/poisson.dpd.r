#' MDPDE for Logistic Regression Model 
#'
#' Computes the robust Minimum Density Power Divergence Estimator (MDPDE) under a Logistic Regression Model 
#'
#' Reference: Ghosh, A., & Basu, A. (2016). Robust estimation in generalized linear models: the density power divergence approach. Test, 25(2), 269-290.
#'
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
#' @param alpha Numeric. The DPD tuning parameter (0<= alpha <=1)
#' @param Method String. Numerical optimization method to be used for the estimation process. Possible options are "L-BFGS-B","Nelder-Mead", "BFGS", "CG", which are the same as the input of 'optim' function in R. Optional. Default is "L-BFGS-B".
#' @param p Integer. Number of columns in the covariate matrix (X). Optional.
#' @param Initial Vector. Initial values of the parameters for the estimation process (initial of sigma needs to be given via log(sigma)). Optional. Default is 1 for all slope parameters, median of y for intercept and MAD(y) for sigma.

#' @return A List of two elements.
#' The first list element is a [p X 1] vector containing the parameter estimates.
#' The second list element is an indicator if the numerical optimization within the estimation process has converged: 0 == convergence. (Same as the convergence from 'optim' function of R)
#'
#'
#' @examples
#' \donttest{
#' n <- 50; p <- 5;
#' beta <- rep(1,p); sigma<-1
#' SigmaX <- diag(p-1);
#' X <- mvrnorm(n, mu=rep(0,p-1), Sigma=SigmaX);
#' X0 <- cbind(1,X)
#' Y <- drop(X0 %*% beta + sigma*rnorm(n))
#'
#' alpha <- 0.3
#' MDPDE<-lmdpd(Y,X0,alpha)
#' }
#'
#' @export


poisson.dpd <- function(y,X,alpha,method = "L-BFGS-B",p=ncol(X),Initial=matrix(rep(1,p))) {

  n <- length(y)
  T<-solve(t(X)%*%X)

#Objective function for the computation of MDPDE (Minimum DPD estimators)
objfunction_Poiss_reg <- function(t) {
    fa<-matrix(data=NA,nrow=n,ncol=1)
    fb<-matrix(data=NA,nrow=n,ncol=1)
    fc<-matrix(data=NA,nrow=n,ncol=1)
    
    
    for (i in 1:n) {
      lambda<-exp(X[i,]%*%t)
      fa[i,] <- lambda - (y[i] *(X[i,] %*%t))
      fb[i,] <- dpois(y[i],lambda)^alpha
      s <- rpois(10000,lambda)
      fc[i,]=mean((dpois(s,lambda)^alpha))
    }
    if ( alpha == 0) 
    { f = sum(fa)}
    else
    { f =sum(fc)-((1+alpha)*sum(fb)/alpha)}
  }
  


## Obtain the MDPDE in 'est'
## Default method used in optim is "L-BFGS-B", in case of non-convergence, this method can be changed
  result<-optim(Initial, objfunction_Poiss_reg , gr = NULL,
              method =  method,
              lower = -Inf, upper = Inf, control = list(), hessian = FALSE)
  est=result$par

## Record convergence indicator in 'conv'
  conv=result$convergence

## Output
  output<-list(est, conv)
  return(output)
}
