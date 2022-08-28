#' MDPDE for Linear Regression Model with knwon error variance (equal to one)
#'
#' Computes the robust Minimum Density Power Divergence Estimator (MDPDE) under a Linear Regression Model  y=X*beta + e, with e ~ N(0, 1).
#'
#' Reference: Ghosh, A., & Basu, A. (2013). Robust estimation for independent non-homogeneous observations using density power divergence with applications to linear regression. Electronic Journal of statistics, 7, 2420-2456.
#'
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
#' @param alpha Numeric. The DPD tuning parameter (0<= alpha <=1)
#' @param Method String. Numerical optimization method to be used for the estimation process. Possible options are "L-BFGS-B","Nelder-Mead", "BFGS", "CG", which are the same as the input of 'optim' function in R. Optional. Default is "L-BFGS-B".
#' @param p Integer. Number of columns in the covariate matrix (X). Optional.
#' @param Initial Vector. Initial values of the parameters for the estimation process (initial of sigma needs to be given via log(sigma)). Optional. Default is 1 for all slope parameters, median of y for intercept and MAD(y) for sigma.

#' @return A List of two elements.
#' The first list element is a [(p+1) X 1] vector containing the parameter estimates (last entry being estimate of sigma)
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


lrm.dpd <- function(y,X,alpha,method = "L-BFGS-B",p=ncol(X),Initial=matrix(rep(1,p))) {

  n <- length(y)
  T<-solve(t(X)%*%X)

#Objective function for the computation of MDPDE (Minimum DPD estimators)
objfunction_Linear_reg <- function(t) {
    beta<-t[1:p]
    sigma<-1
    r<- (-((y-X%*%beta)^2)/(2*sigma*sigma))
    c<-1/(sigma^alpha)

    if ( alpha == 0)
      { f = -mean(r)+t[p+1]}
    else
      { f =c*((1/sqrt(1+alpha))-((1+alpha)*mean(exp(r))/alpha))}
    return(f)
}


## Obtain the MDPDE in 'est'
## Default method used in optim is "L-BFGS-B", in case of non-convergence, this method can be changed
  result<-optim(Initial, objfunction_Linear_reg , gr = NULL,
              method =  method,
              lower = -Inf, upper = Inf, control = list(), hessian = FALSE)
  est=result$par

## Record convergence indicator in 'conv'
  conv=result$convergence

## Output
  output<-list(est, conv)
  return(output)
}
