# SmoothLP

Some work-in-progress `.ado` files to perform Smooth Local Projection estimation outlined in [Barnichon and Brownlees (2019)](https://www.mitpressjournals.org/doi/abs/10.1162/rest_a_00778), BB19 hereafter. This was written based on the [replication files](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/8KQJBJ) published alongside the paper (available in R and Matlab). A Julia implementation can be found [here](https://github.com/justinjjlee/SmoothLocalProjections.jl). [Li, Plagborg-Møller, and Wolf (2024)](https://www.sciencedirect.com/science/article/pii/S030440762400068X?via%3Dihub) include this estimator in their paper assessesing the performance of variants of local projections and VARs across a litany of simulations ([GitHub repo](https://github.com/dake-li/lp_var_simul)). 

Smooth Local Projections extend penalize regression methods to local projections (LP). The estimator of BB19 comes from minimizing a ridge-type loss function with a penalty matrix constructured from B-spline basis functions. A major motivation is the use of local projections to generate impulse response functions (IRFs). More specifically, the relevant deliverable is a plot of the coeficients on $x_t$ from regressions of $y_{t+h}$ on $(x_t,\boldsymbol{W}_t)$ for $h=h_1...H$, representing the response path of an outcome varible of interest $y$ given a change in $x$ at $t$. When estimated via a default LP structure, these coefficients can have a lot of variability from one horizon to the next (created "jagged" plots), as well as large standard errors. BB19 addresses both, most obviously because nonzero values of the ridge parameter $\lambda$ produce more efficient estimates, but the estimated IRFs shrink to a polynomial as $\lambda$ grows because higher order derivitives with this B-spline penalty matrix converge to 0, making response paths even "smoother" than a more straightforward penalized approach. 

## Replication of Example

Section 3 of the paper shows an application of this Smooth LP estimator through estimates the impulse response of output to interest rates. The codes below reproduce the IRF with 90% confidence bands. 

```
clear all
import delimited "https://raw.githubusercontent.com/paulbousquet/SmoothLP/main/data.csv", clear

* if you want your control set to include lags, set as time series data
gen mdate = q(1959q1) + _n-1
tsset mdate, q

* outcome variable 
local y yg
* endogenous variable of interest 
local x ir
* additional covariates 
local w pi 

* which period the response begins 
local h1 = 1
* horizon of IRF 
local H = 20

* Vector of possible penalization parameters 
local lambda 0.00001 0.0001 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 .009 .01

lproj_irf `y' `x' `w', h(`H') h1(`h1') lambda(`lambda') k(5) lag(4)

```
In most Macro setting of interest, the IRF settingg gives rise to an identification issue because of the endogeneity of $x$. One approach, referred to "identification through controls" in BB19, is to assume macro variables of interest evolve according to a VAR system, essentially allowing one to back out an exogenous shock to $x$ by conditioning on covariates $\boldsymbol{W}_t$. As they point out in their paper, for this application with output and interest rates, this can be thought of as identifying shocks in a Taylor Rule. However, this identification through controls approach has grown less popular over time because of the need to achieve identification from stuctural assumptions -- LPs are attactive in the first place because the little structur imposed as a default has shown to yield much less biased results in finite samples than VARs. However, even if we don't feel comfortable ascribing a causaul interpretation to what's being estimated, uncovering correlations is still important and it's quite easy to change the conditioning set to see exactly how sensitive the results are. 

## Syntax 

```
syntax varlist(min=2) [if] [in], H(integer) Lambda(numlist) K(integer) [H1(integer 0)] [R(integer 2)] [Lag(integer 0)]

```
* You can call `mksplit` just like `reg`. To plot the IRF of `y` to `x`, list `y x` in that order. Every variable listed after `x` will be included in the list of controls.
* H is the horizon length
* Lambda is a list of numbers considered for the ridge parameter $\lambda$ 
* K is the number of cross-validation folds
* H1 is the period the IRF starts
* r is the order of the limit polynomial. More specifically, the r+1-th derivitive of the IRF converges to 0 as $\lambda$ grows.
* Lag allows you to include lags of control variables in the conditioning set. For this, a `tsset` command must be run before `lproj_irf`

 
 ## Installation

 Because this is the first time I have tried to create a Stata package (or worked with .ado files in general), for now I am going to hold off on making it available on ssc. Please [email me](mailto:ptb8zf@virginia.edu) if you have any comments or feature requests. I will be updating it this summer as I work through a couple projects. 

 ```
net from https://raw.githubusercontent.com/paulbousquet/SmoothLP/master 
```
