# dpdSIS
Robust Sure Independence Screening (DPD-SIS) and Conditional Screening (DPD-CSIS) using the Minimum Density Power Divergence (DPD) Estimators 
(Ghosh and Thoresen, 2021; Ghosh et al., 2021)


<!-- badges: start -->

<!-- badges: end -->

## Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("abhianik/dpdSIS")
```

## Examples: DPD-SIS and DPD-CSIS for a specified GLM

You first need to load the following packages.

``` r
library(dpdSIS)

## Additional packages Required
library(MASS)
library(parallel)
library(doParallel)
library(sfsmisc)

```


Examples to run DPD-SIS given the response vector y and covariate matrix X (without intercept) for a given GLM. 

``` r
n = length(y)          # No. of observations
d = floor(n/log(n))     # Required model size
alpha = 0.1             # DPD tuning parameter

## Run DPD-SIS Screening for a given alpha and specified GLM

SIS <- dpd.sis(d, y, X, alpha, reg='lrm')        # For Linear regression model (with known error variance equal to one)
SIS <- dpd.sis(d, y, X, alpha, reg='logistic')   # For Logistic regression model
SIS <- dpd.sis(d, y, X, alpha, reg='poisson')    # For Poisson regression model

SIS <- dpd.sis.lrm(d, y, cbind(1, X), alpha)        # For Linear regression model with unknwon error variance

```


Examples to run DPD-CSIS given the response vector y and covariate matrix X (without intercept) for a given GLM, when a matrix of conditioning covarite is given as XC (shoudl include intercept) 

``` r
n = length(y)           # No. of observations
d = floor(n/log(n))     # Required model size
alpha = 0.1             # DPD tuning parameter

## Run DPD-SIS Screening for a given alpha and specified GLM

SIS <- dpd.sis(d, y, X, alpha, reg='lrm', XC)        # For Linear regression model
SIS <- dpd.sis(d, y, X, alpha, reg='logistic', XC)   # For Logistic regression model
SIS <- dpd.sis(d, y, X, alpha, reg='poisson', XC)    # For Poisson regression model

```

The function dpd.sis (and also dpd.sis.lrm) returns a matrix (SIS) of order d X 2, where the first column contains the ranked variable indeces (in decreasing order of abosolute values of the marginal slope) and the second column contains the estimated marginal slopes of the corresponding variables. It returns the top d many variables in decreasing order of importance. 

## The usual SIS methods (Fan and Lv, 2008)

``` r
library(SIS)

X0 = cbind(1, X)
SIS_usual <- SIS(X0, y, family='gaussian', iter = FALSE, nsis=d)

# Indices of screened variables are to be obtained as follows: 
I1=SIS_usual$sis.ix0

```


## Competitive existing robust SIS methods, available ONLY for Linear Regression Model!

In the following, the covarite matrix X0 should include the intercept variable as well, i.e., X0 = cbind(1, X), and the output obtained in each case has the exact same structure as the output of dpd.sis or dpd.sis.lrm (for ease of comparison).

  - Rank-SIS: (Li et al., 2012)

<!-- end list -->

``` r
SIS_rank <- rank.sis.lrm(d, y, X0)

```

  - GK-SIS: (Gather and Guddat, 2008)

<!-- end list -->

``` r
SIS_GK <- GK.sis.lrm(d, y, X0)

```

  - dCorr-SIS: (Wan et al., 2017)

<!-- end list -->

``` r
library(energy)

SIS_dCor <- dcor.sis.lrm(d, y, X0)

```

  - MCP-SIS: (Mu and Xiong, 2014)

<!-- end list -->

``` r
SIS_MCP <- MCP.sis.lrm(d, y, X0)

```


## References 

Fan J, Lv J. Sure independence screening for ultrahigh dimensional feature space. J Royal Stat Soc B. 2008; 70(5):849{911.

Gather U, Guddat C. Comment on "Sure Independence Screening for Ultrahigh Dimensional Feature Space" by Fan, JQ and Lv, J. J Royal Stat Soc B. 2008; 70:893-895.

Ghosh A, Thoresen M. A robust variable screening procedure for ultra-high dimensional data. Stat Meth Med Res, 2021, 30(8), 1816–1832.

Ghosh A, Ponzi E, Sandanger T,  Thoresen M. Robust Sure Independence Screening for Non-polynomial dimensional Generalized Linear Models. arXiv preprint 2021; arXiv:2005.12068v2.

Li G, Peng H, Zhang J, Zhu L. Robust rank correlation based screening. Ann Stat. 2012; 40(3):1846-1877.

Mu W, Xiong S. Some notes on robust sure independence screening. J App Stat. 2014; 41(10):2092-2102.

Wang T, Zheng L, Li Z, Liu H. A robust variable screening method for high-dimensional data. J App Stat. 2017; 44(10):1839-1855.
