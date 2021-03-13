## R functions for computing several existing non-robust and robust variable screening for linear regression model (LRM)



#' Rank correlation based SIS under a Linear Regression Model
#'
#' Performs the rank correlation based Robust Variable Screening (DPD-SIS) under a Linear Regression Model  y=X*beta + e, with e ~ N(0, sigma^2) for some unknown sigma, using parallel computation.
#'
#' Reference: Li G, Peng H, Zhang J,  Zhu L.   Robust rank correlation based screening. {Ann Stat}. 2012; 40(3):1846-1877.
#'
#' @param d Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
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
#' SIS<-rank.sis.lrm(d,Y,X0)
#' }
#'
#' @export


rank.sis.lrm <- function(d,y,X) {

  n <- length(y)
  p=ncol(X)
  # T<-solve(t(X)%*%X)
  #Output matrices
  est<-matrix(0,p-1,1);


 # for(w in 1:(p-1)){
library("doParallel")
 nprocs=detectCores() - 1
 cl=makeCluster(nprocs)
 registerDoParallel(cl)
   est<-foreach(w=1:(p-1),.combine=rbind,.packages = c("sfsmisc","MASS")) %dopar% {

  y=as.numeric(y)
  X1=X[,(w+1)]

  est[w]=cor(y,X1,method="spearman")
  #output for parallel processing Out<-est
}#close for cycle w
   stopCluster(cl)

I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
I1=I$ix[1:d]
beta=est[I1]
result=cbind(I1, beta)
return(result)
}


###########################################################################################################################


#' Gnanadesikan-Kettenring correlation based SIS under a Linear Regression Model
#'
#' Performs the Gnanadesikan-Kettenring correlation based Robust Variable Screening (DPD-SIS) under a Linear Regression Model  y=X*beta + e, with e ~ N(0, sigma^2) for some unknown sigma, using parallel computation.
#'
#' Reference: Gather U, Guddat C. Comment on ``Sure Independence Screening for Ultrahigh Dimensional Feature Space" by Fan, JQ and Lv, J. {J Royal Stat Soc B}. 2008; 70:893-895.
#'
#' @param d Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
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
#' SIS<-GK.sis.lrm(d,Y,X0)
#' }
#'
#' @export

GK.sis.lrm <- function(d,y,X) {

  n <- length(y)
  p=ncol(X)
  # T<-solve(t(X)%*%X)
  #Output matrices
  est<-matrix(0,p-1,1);


  # for(w in 1:(p-1)){
  library("doParallel")
  nprocs=detectCores() - 1
  cl=makeCluster(nprocs)
  registerDoParallel(cl)
  est<-foreach(w=1:(p-1),.combine=rbind,.packages = c("robustbase")) %dopar% {

    y=as.numeric(y)
    X1=X[,(w+1)]

    a=covGK(y,X1)
    est[w]=sqrt(abs(a)/(covGK(y,y)*covGK(X1,X1)))*sign(a)
    #output for parallel processing Out<-est
  }#close for cycle w
  stopCluster(cl)

  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta)
  return(result)
}



###########################################################################################################################


#' Distance correlation based SIS under a Linear Regression Model
#'
#' Performs the distance correlation based Robust Variable Screening (DPD-SIS) under a Linear Regression Model  y=X*beta + e, with e ~ N(0, sigma^2) for some unknown sigma, using parallel computation.
#'
#' Reference: Li R, Zhong W, Zhu L.  Feature screening via distance correlation learning.  {J Amer Statist Assoc}. 2012; 107(499):1129-1139.
#'
#' @param d Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
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
#' SIS<-dcor.sis.lrm(d,Y,X0)
#' }
#'
#' @export

dcor.sis.lrm <- function(d,y,X) {

  n <- length(y)
  p=ncol(X)
  # T<-solve(t(X)%*%X)
  #Output matrices
  est<-matrix(0,p-1,1);


  # for(w in 1:(p-1)){
  library("doParallel")
  nprocs=detectCores() - 1
  cl=makeCluster(nprocs)
  registerDoParallel(cl)
  est<-foreach(w=1:(p-1),.combine=rbind,.packages = c("energy")) %dopar% {

    y=as.numeric(y)
    X1=X[,(w+1)]

    Ys=(y-median(y))/mad(y)
    Xs=(X1-median(X1))/mad(X1)
    a=dcor.test(Ys,Xs)
    est[w]=a$estimate[2]
    #output for parallel processing Out<-est
  }#close for cycle w
  stopCluster(cl)

  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta)
  return(result)
}



###########################################################################################################################


#' MCP correlation based SIS under a Linear Regression Model
#'
#' Performs the MCP correlation based Robust Variable Screening (DPD-SIS) under a Linear Regression Model  y=X*beta + e, with e ~ N(0, sigma^2) for some unknown sigma, using parallel computation.
#'
#' Reference: Mu W, Xiong S.  Some notes on robust sure independence screening.  {J App Stat}. 2014; 41(10):2092--2102.
#'
#' @param d Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)
#' @param y Vector.   The response vector [n X 1].
#' @param X Matrix. Covariate Matrix  [n X p].
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
#' SIS<-MCP.sis.lrm(d,Y,X0)
#' }
#'
#' @export

MCP.sis.lrm <- function(d,y,X) {

  n <- length(y)
  p=ncol(X)
  # T<-solve(t(X)%*%X)
  #Output matrices
  est<-matrix(0,p-1,1);


  # for(w in 1:(p-1)){
  library("doParallel")
  nprocs=detectCores() - 1
  cl=makeCluster(nprocs)
  registerDoParallel(cl)
  est<-foreach(w=1:(p-1),.combine=rbind,.packages = c("MASS")) %dopar% {

    y=as.numeric(y)
    X1=X[,(w+1)]

    Ys=(y-median(y))/mad(y)
    Xs=(X1-median(X1))/mad(X1)
    z=Ys*Xs
    est[w]=median(z)
    #output for parallel processing Out<-est
  }#close for cycle w
  stopCluster(cl)

  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta)
  return(result)
}
