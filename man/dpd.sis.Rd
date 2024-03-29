% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpd.sis.r
\name{dpd.sis}
\alias{dpd.sis}
\title{DPD-SIS and DPD-CSIS under a Generalized Linear Model (GLM) with no variance parameter}
\usage{
dpd.sis(d, y, X,  alpha,  reg = c("lrm", "logistic", "poisson"),  XC = matrix(rep(1, nrow(X))), Initial = matrix(rep(1, ncol(X))), Method = "L-BFGS-B")
}
\arguments{
\item{d}{Integer. A list of the data matrices. Number of features to be selected. (1<= d <=p)}

\item{y}{Vector.   The response vector [n X 1].}

\item{X}{Matrix. Covariate Matrix  [n X p]. It shoudl only include variables to be considered for screening. So, it should not contain the intercept or any conditioning variables.}

\item{alpha}{Numeric. The DPD tuning parameter (0<= alpha <=1)}

\item{reg}{A string indiacting the regression model. Possible options are  "lrm" (default), "logistic"  and "poisson", indicating the linear regression (with unit error variance), logistic regression and poisson regression, respectively.}

\item{Initial}{Vector. Initial values of the marginal slope parameter for estimation process. Default: 1 for all parameters. Optional.}

\item{Method}{String. Numerical optimization method to be used for computation of marginal slopes. Possible options are "L-BFGS-B","Nelder-Mead", "BFGS", "CG", which are the same as the input of 'optim' function in R. Optional. Default is "L-BFGS-B".}

\item{X_C}{Matrix. Conditioning Covariate Matrix [n X q] in DPD-CSIS. Default: only the intercept variable (q=1). Optional for DPD-SIS.}
}
\value{
SIS Vector. A d X 2 vector where the first column contains the variable index (ranked in decreasing order of importance) and the second column consist of the corresponding MDPDE of the slope
}
\description{
Performs the Density Power Divergence based Robust Variable Screening (DPD-SIS) and associated Robust Conditional Variable Screening (DPD-CSIS) under a GLM with no variance parameter, using parallel computation.
}
\details{
Reference: Ghosh A, Ponzi E, Sandanger T, Thoresen M. Robust Sure Independence Screening for Non-polynomial dimensional Generalized Linear Models. arXiv preprint 2021; arXiv:2005.12068v2.
}
\examples{
\donttest{
n <- 50; p <- 500;
beta <- rep(0,p); beta[c(1:5)] <- c(1,1,1,1,1);
d <- floor(n/log(n));
Sigma <- diag(p-1);
X <- mvrnorm(n, mu=rep(0,p-1), Sigma=Sigma);
X0 <- cbind(1,X)
Y <- drop(X0 \%*\% beta + 2*rnorm(n))

alpha <- 0.3
SIS<-dpd.sis(d,Y,X,alpha, reg='lrm')
}

}
