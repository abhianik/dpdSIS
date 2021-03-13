# dpdSIS
Robust Sure Independence Screening (SIS) using the Minimum Density Power Divergence (DPD) Estimators 
(Ghosh and Thoresen, 2020; Ghosh et al., 2021)


<!-- badges: start -->

<!-- badges: end -->

## Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("abhianik/dpdSIS")
```

## Example

An example to run DPD-SIS for a simulated dataset from linear regression model:


``` r
library(dpdSIS)

## Additional Library Required
library(MASS)
library(parallel)
library(doParallel)
library(sfsmisc)

## Simulate the data
n=50;     # No. of observations
p=5000;   # No. of variables
beta = rep(0,p)
beta[c(1:5)] = c(1,1,1,1,1); # True non-zero coefficients
  
d=floor(n/log(n))     # Required model size
set.seed(rep1*1e2)

Sigma = diag(p-1)
#Sigma[Sigma==0]=0.5
X = mvrnorm(n, mu=rep(0,p-1), Sigma=Sigma) # Simulate covariates
X0 = cbind(1,X)

Y = drop(X0 %*% beta + 2*rnorm(n))  # Simulate reponse
    
## add outliers in response
epsilon=0.1                   # outlier proportion
nout = ceiling(n*epsilon)     # outlier numbers
Y[1:nout] = Y[1:nout] - 30    # adding outliers

    
## Run DPD-SIS Screening for a given alpha 
  
alpha=0.3
SIS<-dpd.sis.lrm(d,Y,X0,alpha) 
```

The function returns a matrix (SIS) of order d X 2, where the first column contains the ranked variable indeces (in decreasing order of abosolute values of the marginal slope) 
and the second column contains the estimated marginal slopes of the corresponding variables. It returns the top d many variables in decreasing order of importance. 


The following codes illustrate the compuation of different summary measures from the output (with d=p)

``` r
I1=SIS0[(1:(n-1)),1]
I2=SIS0[,1]

# True positives
TP=sum(I1==1)+sum(I1==2)+sum(I1==3)+sum(I1==4)   

# Indicator if the full model is correctly selected
full=as.numeric(tp==4)                           

# Minimum model size (MMS) required to cover the correct model fully
minM=max(which(I2==1), which(I2==2), which(I2==3), which(I2==4))  
```

## Comparison with the SIS and other existing robust SIS methods

  - Usual SIS: (Fan and Lv, 2008)

<!-- end list -->

``` r
library(SIS)

SIS_usual<- SIS(X0, Y, family='gaussian', iter = FALSE, nsis=n-1)

I1=SIS_usual$sis.ix0
TP=sum(I1==1)+sum(I1==2)+sum(I1==3)+sum(I1==4)
full=as.numeric(tp==4)
    
    
```

  - Rank-SIS: (Li et al., 2012)

<!-- end list -->

``` r
SIS_rank<-rank.sis.lrm(p,Y,X0)

I1=SIS_rank[(1:(n-1)),1]
I2=SIS_rank[,1]
TP=sum(I1==1)+sum(I1==2)+sum(I1==3)+sum(I1==4)   
full=as.numeric(tp==4)                           
minM=max(which(I2==1), which(I2==2), which(I2==3), which(I2==4))  
```

  - GK-SIS: (Gather and Guddat, 2008)

<!-- end list -->

``` r
SIS_GK<-GK.sis.lrm(p,Y,X0)

I1=SIS_GK[(1:(n-1)),1]
I2=SIS_GK[,1]
TP=sum(I1==1)+sum(I1==2)+sum(I1==3)+sum(I1==4)   
full=as.numeric(tp==4)                           
minM=max(which(I2==1), which(I2==2), which(I2==3), which(I2==4))  
```

  - dCorr-SIS: (Wan et al., 2017)

<!-- end list -->

``` r
library(energy)

SIS_dCor<-dcor.sis.lrm(p,Y,X0)

I1=SIS_dCor[(1:(n-1)),1]
I2=SIS_dCor[,1]
TP=sum(I1==1)+sum(I1==2)+sum(I1==3)+sum(I1==4)   
full=as.numeric(tp==4)                           
minM=max(which(I2==1), which(I2==2), which(I2==3), which(I2==4))  
```

  - MCP-SIS: (Mu and Xiong, 2014)

<!-- end list -->

``` r
SIS_MCP<-MCP.sis.lrm(p,Y,X0)

I1=SIS_MCP[(1:(n-1)),1]
I2=SIS_MCP[,1]
TP=sum(I1==1)+sum(I1==2)+sum(I1==3)+sum(I1==4)   
full=as.numeric(tp==4)                           
minM=max(which(I2==1), which(I2==2), which(I2==3), which(I2==4))  
```


## References 

Fan J, Lv J. Sure independence screening for ultrahigh dimensional feature space. J Royal Stat Soc B. 2008; 70(5):849{911.

Gather U, Guddat C. Comment on "Sure Independence Screening for Ultrahigh Dimensional Feature Space" by Fan, JQ and Lv, J. J Royal Stat Soc B. 2008; 70:893-895.

Ghosh A, Thoresen M. A robust variable screening procedure for ultra-high dimensional data. arXiv preprint. 2021; arXiv:2004.14851.

Ghosh A, Ponzi E, Sandanger T,  Thoresen M. Robust Sure Independence Screening for Non-polynomial dimensional Generalized Linear Models. arXiv preprint 2021; arXiv:2005.12068v2.

Li G, Peng H, Zhang J, Zhu L. Robust rank correlation based screening. Ann Stat. 2012; 40(3):1846-1877.

Mu W, Xiong S. Some notes on robust sure independence screening. J App Stat. 2014; 41(10):2092-2102.

Wang T, Zheng L, Li Z, Liu H. A robust variable screening method for high-dimensional data. J App Stat. 2017; 44(10):1839-1855.
