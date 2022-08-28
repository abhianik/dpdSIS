#' DPD-SIS and DPD-CSIS under a Generalized Linear Model (GLM) with no variance parameter 
#'
#' Performs the Density Power Divergence based Robust Variable Screening (DPD-SIS) and associated Robust Conditional Variable Screening (DPD-CSIS) under a GLM with no variance parameter, using parallel computation.
#'
#' Reference: Ghosh A, Ponzi E, Sandanger T, Thoresen M. Robust Sure Independence Screening for Non-polynomial dimensional Generalized Linear Models. arXiv preprint 2021; arXiv:2005.12068v2.
#'
#' @param d Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p]. It shoudl only include variables to be considered for screening. So, it should not contain the intercept or any conditioning variables. 
#' @param alpha Numeric. The DPD tuning parameter (0<= alpha <=1)
#' @param reg. A string indiacting the regression model. Possible options are  "lrm" (default), "logistic"  and "Poisson", indicating the linear regression (with unit error variance), logistic regression and poisson regression, respectively. 
#' @param X_C Matrix. Conditioning Covariate Matrix [n X q] in DPD-CSIS. Default: only the intercept variable (q=1). Optional for DPD-SIS.
#' @param Initial Vector. Initial values of the marginal slope parameter for estimation process. Default: 1 for all parameters. Optional.
#' @param Method String. Numerical optimization method to be used for computation of marginal slopes. Possible options are "L-BFGS-B","Nelder-Mead", "BFGS", "CG", which are the same as the input of 'optim' function in R. Optional. Default is "L-BFGS-B".
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
#' SIS<-dpd.sis(d,Y,X,alpha, reg='lrm')
#' }
#'
#' @export

dpd.sis <- function(d, y, X, alpha, reg=c('lrm', 'logistic', 'lrm', 'poisson'), XC=matrix(rep(1,nrow(X))), Initial=matrix(rep(1,ncol(X))), Method="L-BFGS-B") {
  
  n <- length(y)
  p=ncol(X)
  p1=ncol(XC)
  f=match.fun(paste(reg,'.dpd',sep=""))
  est<-matrix(0,p,1);


#### fitting the ROBUST DPD based GLM for given alpha -------

library("doParallel")
nprocs=detectCores();
if (nprocs>1) nprocs=nprocs-1;
cl=makeCluster(nprocs)
registerDoParallel(cl)


est<-foreach(w=1:p,.combine=rbind,.packages = c("sfsmisc","MASS")) %dopar% {
  source(paste(reg,'.dpd.r',sep=""));
  y=as.numeric(y)
  X1=cbind(XC, X[,w])
  init=matrix(c(rep(1,ncol(XC)),Initial[w]))
  
  ldpd<-list(NA*init,0)
  tryCatch(
  ldpd<-f(y,X1,alpha,Method, Initial=init),
  error=function(e){print(e)}
  )
  paste(ldpd[[2]])
  betaw=ldpd[[1]]
  est[w]=betaw[p1+1]
  #output for parallel processing Out<-est
}#close for cycle w

stopCluster(cl)

I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
I1=I$ix[1:d]
beta=est[I1]
result=cbind(I1, beta)
return(result)
}
