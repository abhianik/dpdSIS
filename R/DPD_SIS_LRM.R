#' DPD-SIS under a Linear Regression Model
#'
#' Performs the Density Power Divergence based Robust Variable Screening (DPD-SIS) under a Linear Regression Model  y=X*beta + e, with e ~ N(0, sigma^2) for some unknown sigma, using parallel computation.
#'
#' Reference: Ghosh A, Thoresen M. A robust variable screening procedure for ultra-high dimensional data. arXiv preprint. 2021; arXiv:2004.14851.
#'
#' @param d Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
#' @param alpha Numeric. The DPD tuning parameter (0<= alpha <=1)
#' @param Method String. Numerical optimization method to be used for computation of marginal slopes. Possible options are "L-BFGS-B","Nelder-Mead", "BFGS", "CG", which are the same as the input of 'optim' function in R. Optional. Default is "L-BFGS-B".
#' @param p Integer. Number of columns in the covariate matrix (X). Optional.
#' @param Initial Vector. Initial values of the marginal slope parameter for estimation process. Default: 1 for all parameters. Optional.
#'
#' @return SIS Vector. A d X 2 vector where the first column contains the variable index (ranked in decreasing order of importance) and the second column consist of the corresponding MDPDE of the slope
#'
#' @examples
#' \donttest{
#' n <- 50; p <- 500;
#' beta <- rep(0,p); beta[c(1:5)] <- c(1,1,1,1,1);
#' d <- floor(n/log(n));
#' Sigma <- diag(p-1);
#' X <- mvrnorm(n, mu=rep(0,p-1), Sigma=Sigma);
#' X0 <- cbind(1,X)
#' Y <- drop(X0 %*% beta + 2*rnorm(n))
#'
#' alpha <- 0.3
#' SIS<-dpd.sis.lrm(d,Y,X0,alpha)
#' }
#'
#' @export

dpd.sis.lrm <- function(d,y,X,alpha,Method="L-BFGS-B",p=ncol(X),Initial=matrix(rep(1,p))){

  n <- length(y)
  p=ncol(X)
  # T<-solve(t(X)%*%X)
  est<-matrix(0,p-1,1);


#### fitting the ROBUST DPD based LINEAR MODEL for given alpha -------

# Initial<-matrix(c(mu, rep(0,ncol(X)-1),0.5)) #Initial value for computing DPD estimate
# for(w in 1:(p-1)){
library("doParallel")
nprocs=detectCores() - 2
cl=makeCluster(nprocs)
registerDoParallel(cl)
est<-foreach(w=1:(p-1),.combine=rbind,.packages = c("sfsmisc","MASS")) %dopar% {


  y=as.numeric(y)
  X1=X[,c(1,w+1)]
  init=matrix(c(median(y), Initial[w+1],log(mad(y))))

  ldpd<-lmdpd(y,X1,alpha,Method,2,init)
  paste(ldpd[[2]])
  betaw=ldpd[[1]]
  est[w]=betaw[2]
  #output for parallel processing Out<-est
}#close for cycle w

stopCluster(cl)

I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
I1=I$ix[1:d]
beta=est[I1]

result=cbind(I1, beta)
return(result)
}
