% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic.dpd.r
\name{logistic.dpd}
\alias{logistic.dpd}
\title{MDPDE for Logistic Regression Model}
\usage{
logistic.dpd(y, X, alpha, method = "L-BFGS-B", p = ncol(X), Initial = matrix(rep(1, p)))
}
\arguments{
\item{y}{Vector.   The response vector [n X 1].}

\item{X}{Matrix. Covariate Matrix  [n X p].}

\item{alpha}{Numeric. The DPD tuning parameter (0<= alpha <=1)}

\item{p}{Integer. Number of columns in the covariate matrix (X). Optional.}

\item{Initial}{Vector. Initial values of the parameters for the estimation process (initial of sigma needs to be given via log(sigma)). Optional. Default is 1 for all slope parameters, median of y for intercept and MAD(y) for sigma.}

\item{Method}{String. Numerical optimization method to be used for the estimation process. Possible options are "L-BFGS-B","Nelder-Mead", "BFGS", "CG", which are the same as the input of 'optim' function in R. Optional. Default is "L-BFGS-B".}
}
\value{
A List of two elements.
The first list element is a [p X 1] vector containing the parameter estimates.
The second list element is an indicator if the numerical optimization within the estimation process has converged: 0 == convergence. (Same as the convergence from 'optim' function of R)
}
\description{
Computes the robust Minimum Density Power Divergence Estimator (MDPDE) under a Logistic Regression Model
}
\details{
Reference: Ghosh, A., & Basu, A. (2016). Robust estimation in generalized linear models: the density power divergence approach. Test, 25(2), 269-290.
}
\examples{
\donttest{
n <- 50; p <- 5;
beta <- rep(1,p); sigma<-1
SigmaX <- diag(p-1);
X <- mvrnorm(n, mu=rep(0,p-1), Sigma=SigmaX);
X0 <- cbind(1,X)
pr = 1/(1+exp(-(X0 \%*\% beta)))
Y = rbinom(n,1,pr)

alpha <- 0.3
MDPDE<-logistic.dpd(Y,X0,alpha)
}

}
