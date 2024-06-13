# SmoothLP

Some work-in-progress `.ado` files for Smooth Local Projection estimation outlined in [Barnichon and Brownlees (2019)](https://www.mitpressjournals.org/doi/abs/10.1162/rest_a_00778), BB19 hereafter. This was written based on their [replication files](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/8KQJBJ) (available in R and Matlab). A Julia implementation can be found [here](https://github.com/justinjjlee/SmoothLocalProjections.jl). [Li, Plagborg-Møller, and Wolf (2024)](https://www.sciencedirect.com/science/article/pii/S030440762400068X?via%3Dihub) assess the performance of of variants of VARs and local projections, including this estimator, across a litany of simulations ([GitHub repo](https://github.com/dake-li/lp_var_simul)). 



<p align="center">
  <a href="#replication-of-example">Replication</a> |
  <a href="#iv-extension">IV</a> |
  <a href="#syntax">Syntax</a> |
  <a href="#future-development">Development</a> |
  <a href="#installation">Installation</a>
</p>



***

Smooth Local Projections extend penalized regression methods to local projections (LP). The estimator of BB19 comes from minimizing a ridge-type loss function with a penalty matrix constructured from B-spline basis functions. A major motivation is the use of local projections to generate impulse response functions (IRFs). More specifically, the relevant deliverable is a plot of the coeficients on $x_t$ from regressions of $y_{t+h}$ on $(x_t,\boldsymbol{W}_t)$ for $h=h_1,\dots,H$, representing the response path of an outcome varible of interest $y$ given a change in $x$ at $t$. When estimated via a default LP structure, these coefficients can have a lot of variability from one horizon to the next (creating "jagged" plots), as well as large standard errors. BB19 addresses both, most obviously because nonzero values of the ridge parameter $\lambda$ produce more efficient estimates, but the estimated IRFs shrink to a polynomial as $\lambda$ grows because higher order derivitives with this B-spline penalty matrix converge to 0, making response paths even "smoother" than a more straightforward penalized approach. 

## Replication of Example

Section 3 of the paper shows an application of estimating the impulse response of output to interest rates. The codes below reproduce the IRF with 90% confidence bands, the primary graph from the `demo` code in their replication files. 

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

* set vmat option to nw for Newey-West, Huber-White by default 

lproj_irf `y' `x' `w', h(`H') h1(`h1') lambda(`lambda') k(5) lag(4) vmat("nw")

```
If you want to customize the graphs, the IRF values (`results1`) and bands (`irc1`, `irc2`) are stored as variables (with `time` as the x axis). For instance, this is what's run as a default, but you can run it on its own after executing `lproj_irf`, as shown below. I should also note that these values are not identical the output from the R replication files, but the difference is subtle enough that I'm tenatively assuming it comes down to well-known precision issues with Mata (e.g., 2/3 has an error after the 7th digit) rather than a transcription error, but please parse for yourself! 
```
tw (rarea irc1 irc2 time, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter result1 time, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if time<=`H'
```

In most Macro settings, this IRF framework gives rise to an identification issue because of the endogeneity of $x$. One approach, referred to "identification through controls" in BB19, is to assume macro variables of interest evolve according to a VAR system, essentially allowing one to back out an exogenous shock to $x$ by conditioning on covariates $\boldsymbol{W}_t$. As they point out, with respect to the application with output and interest rates, this can be thought of as identifying shocks in a Taylor Rule. However, this identification through controls approach has grown less popular over time because of the reliance on stuctural assumptions -- LPs are attactive in the first place in part because the minimal structure imposed has shown to yield much less biased results in finite samples compared to VARs. However, even if we don't feel comfortable ascribing a causaul interpretation to what's being estimated, uncovering correlations is still important and it's quite easy to change the conditioning set to get an idea of how sensitive the results are. 

## IV Extension 

Section 4 of BB19 gives an illustration of the estimator in concert with a LP-IV framework. While no replication files exists for this exercise, conceptually the implementation is straightforward: 2SLS but the first stage contains no regularization procedure. In other words, our `x` becomes the fitted values from the first stage regression. 

To do IV estimation, simply include the endogenous and exogenous variables as you would using `ivregress`. For instance, recall the previous example. If we instrument `ir` with `rr` 

```
lproj_irf `y' `w' (ir=rr), h(`H') h1(`h1') lambda(`lambda') k(5) lag(4) vmat("nw")
```
The code is deisgned to handle multiple instruments. If you want graphs for all endogenous variables that are instrumented, add the `mult` option  

## Syntax 

```
syntax anything(equalok) [if] [in], H(integer) Lambda(numlist) K(integer) [H1(integer 0)] [R(integer 2)] [Lag(integer 0) vmat(string) se(integer 1) CUM MULT NODRAW NOADJ]
```
* You can call `lproj` just like `reg`. To plot the IRF of `y` to `x`, list `y x` in that order. Every variable listed after `x` will be included in the list of controls. See section above for IV option
* H is the horizon length
* Lambda is a list of numbers considered for the ridge parameter $\lambda$ 
* K is the number of cross-validation folds
* H1 is the period the IRF starts
* r sets order of the limit polynomial. More specifically, the $r$-th derivitive of the IRF converges to 0 as $\lambda$ grows 
* Lag allows you to include lags of control variables in the conditioning set. For this, a `tsset` command must be run before `lproj_irf`
* vmat takes option "nw" if you would rather use Newey-West standard errors over Huber-White, but note that [Herbst and Johannsen (2024)](http://www.sciencedirect.com/science/article/pii/S0304407624000010) finds the NW variance matrix will often be biased while [Plagborg-Møller and Montiel Olea (2021)](https://joseluismontielolea.com/lp_inference_ecta.pdf) show that if a sufficient number of lags are included as controls, the usual HW errors are unbiased and autocorrelation robust.
* se allows you to scale the size of the shock (by default, it's a 1 std shock).
* Alternativly, you can add the `noadj` option for a pure plot of the coefficients
* Add `cum` to instead do cumulative IRFs
* If you have multiple IVs and you want IRFs for all endogenous variables, add `mult`
* Add `nodraw` for no graphs

To recap and list all of the things that are stored following the command 
* matrix results: the coeficients that are plotted
* matrix ic: the bands on the coeficients
* results and ic are also stored as variables, along with lagged variables if that option is chosen 
* matrix se: the standard errors
* scalar min_lambda: the $\lambda$ parameter chosen through cross-validation 

## Future Development 

Please [email me](mailto:ptb8zf@virginia.edu) if you have any comments or feature requests. I will be updating it this summer as I work through a couple projects. Some things I have planned 

* more customizability: confidence bands, graphs, ability to have different lag lengths across controls  
* error messages that are common practice in stata packages (e.g., if Lag() is specified but data not loaded with time series format) 
 
 ## Installation

Because this is the first time I have tried to create a Stata package (or worked with .ado files in general), I am going to wait to make it available on ssc while I keep working on it and get feedback. 

There are two dependencies that must be installed first 

 ```
ssc install moremata bspline
net from https://raw.githubusercontent.com/paulbousquet/SmoothLP/master 
```

***
Thank you to my RA, Claud, for all their help!

