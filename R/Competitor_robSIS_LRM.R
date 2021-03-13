#############
## R functions for computing several existing non-robust and robust variable screening for linear regression model (LRM)


#################
### correlation based usual SIS  
## Ref: Fan J, Lv J.  Sure independence screening for ultrahigh dimensional feature space.  {J Royal Stat Soc B }. 2008; 70(5):849--911.


sis.lrm <- function(d,y,X) {
  
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
    
    est[w]=cor(y,X1,method="pearson")
    #output for parallel processing Out<-est
  }#close for cycle w
  stopCluster(cl)
  
  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta)
  return(result)
}



### Rank correlation based SIS
## Ref: Li G, Peng H, Zhang J,  Zhu L.   Robust rank correlation based screening. {Ann Stat}. 2012; 40(3):1846-1877.


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
### gnanadesikan kettenring correlation based SIS 
# Ref: Gather U, Guddat C. Comment on ``Sure Independence Screening for Ultrahigh Dimensional Feature Space" by Fan, JQ and Lv, J. {J Royal Stat Soc B}. 2008; 70:893-895. 

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
### Distance correlation based SIS
## Ref: Li R, Zhong W, Zhu L.  Feature screening via distance correlation learning.  {J Amer Statist Assoc}. 2012; 107(499):1129-1139.


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
### MCP correlation based SIS
## Ref: Mu W, Xiong S.  Some notes on robust sure independence screening.  {J App Stat}. 2014; 41(10):2092--2102.


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
