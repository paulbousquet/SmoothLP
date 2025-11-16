# SmoothLP

Some work-in-progress `.ado` files for Smooth Local Projection estimation outlined in [Barnichon and Brownlees (2019)](https://www.mitpressjournals.org/doi/abs/10.1162/rest_a_00778), BB19 hereafter. This was written based on their [replication files](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/8KQJBJ) (available in R and Matlab). 

* [Li, Plagborg-Møller, and Wolf (2024)](https://www.sciencedirect.com/science/article/pii/S030440762400068X?via%3Dihub) assess the performance of variants of VARs and local projections, including this estimator, across a litany of DGPs ([GitHub repo](https://github.com/dake-li/lp_var_simul))
* A Julia implementation can be found [here](https://github.com/justinjjlee/SmoothLocalProjections.jl) (and [here](https://github.com/junyuan-chen/LocalProjections.jl) as a part of a more comprehesive suite)
* See also an [R package](https://github.com/jackson-mejia/splp/tree/main) for Smooth LPs in a panel data setting 

<p align="center">
  <a href="#installation">Installation</a> |
  <a href="#replication-of-example">Replication</a> |
  <a href="#iv-extension">IV</a> |
  <a href="#inference">Inference </a> |
  <a href="#syntax">Syntax</a> |
  <a href="#future-development">Development</a> 
</p>



***

Smooth Local Projections extend penalized regression methods to local projections (LP). The estimator of BB19 comes from minimizing a ridge-type loss function with a specific penalty matrix motivated by the use of using B-spline basis functions to approximate the coefficient(s) of interest. This approach was conceived to address issues related to a common use of local projections to generate impulse response functions (IRFs). Specifically, the relevant deliverable is a plot of the coeficients on $x_t$ from regressions of $y_{t+h}$ on $(x_t,\bf{W}_t)$ for $h=h_1,\dots,H$, representing the response path of an outcome varible of interest $y$ given a change in $x$ at $t$. When estimated via a default LP structure, these coefficients can have a lot of variability from one horizon to the next (creating "jagged" plots), as well as large standard errors. BB19 addresses both, most obviously because nonzero values of the ridge parameter $\lambda$ produce more efficient estimates, but the estimated IRFs shrink to a polynomial as $\lambda$ grows because higher order derivitives with this penalty matrix in conjunction with the use of B-splines converge to 0, making response paths even "smoother" than a more straightforward penalized approach. Intuitively, smoothness arises across the entire plot because estimating the IRF is treated as a singular optimization problem and the basis functions are constructed to match the number of periods in the IRF (making adjacent coefficients more similar).
 ## Installation
There are two dependencies that must be installed first 

 ```
ssc install moremata; ssc install bspline
copy "https://raw.githubusercontent.com/paulbousquet/SmoothLP/main/slp_irf.ado" slp_irf.ado
```
I am holding off on making this available on ssc while I get more feedback. Please [email me](mailto:pbousquet@virginia.edu) if you have any comments or feature requests (also see the <a href="#future-development"> future development</a>  section).

## Replication of Example

Section 3 of the paper shows an application of estimating the impulse response of output to interest rates. The codes below reproduce the IRF, the primary graph from the `demo` code in their replication files. 

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
* In the original R code, 80% bands are used. Default is 90% 
local ztail = .1
* Inference not well defined in light of cross-validation.
* The example does not account for this and uses a "naive" (unadjusted) starndard error.  
* The default behavior undersmoothes, discussed more in "Inference" section of README
local naive = 1 

* Vector of possible penalization parameters 
local lambda 0.00001 0.0001 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 .009 .01

* set vmat option to nw for Newey-West, Huber-White by default
* discussed in "Inference" section of README

slp_irf `y' `x' `w', h(`H') h1(`h1') lambda(`lambda') k(5) lag(4) vmat("nw") usmooth(`naive') ztail(`ztail')

```
If you want to customize the graphs, the IRF values (`results1`) and bands (`irc1`, `irc2`) are stored as variables (with `time` as the x axis). For instance, this is what's run as a default, but you can run it on its own after executing `slp_irf`, as shown below. 
```
tw (rarea irc1 irc2 time, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter result1 time, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 255))) if time<=`H'
```

In many settings, this specification can be problematic because of the endogeneity of $x$. One approach, referred to "identification through controls" in BB19, is to assume macro variables of interest evolve according to a VAR system, essentially allowing one to back out an exogenous shock to $x$ by conditioning on covariates $\bf{W}_t$. As they point out, with respect to the application with output and interest rates, this can be thought of as identifying shocks in a Taylor Rule. However, this identification through controls approach has grown less popular over time because of the reliance on stuctural assumptions -- LPs are attactive to begin with in part because the minimal structure imposed has been shown to yield much less biased results in finite samples compared to VARs.[^1][^2] But even if we don't feel comfortable ascribing a causal interpretation to what's being estimated, uncovering correlations is still important and it's quite easy to change the conditioning set to get an idea of how sensitive the results are. 

[^1]: For a great practitioner's guide for deciding between LP vs. VAR, see this recent paper from [Montiel Olea, Qian,Plagborg-Møller, and Wolf](https://arxiv.org/abs/2503.17144). They give LPs the edge in general, in large part because any efficiency gains are mechanically a result of insufficient coverage.  
[^2]: [Kolesár and Plagborg-Møller (2025)](https://www.mikkelpm.com/files/nonlinear_causal.pdf) show identification through controls is sensitive to model mispecification/nonlinearity. 
## IV Extension 

Section 4 of BB19 gives an illustration of the estimator in concert with a LP-IV framework. While no replication files exists for this exercise, conceptually the implementation is straightforward: 2SLS but the first stage contains no regularization procedure. In other words, our `x` becomes the fitted values from the first stage regression. 

To do IV estimation, simply include the endogenous and exogenous variables as you would using `ivregress`. For instance, recall the previous example. If we instrument `ir` with `rr` 

```
slp_irf `y' `w' (ir=rr), h(`H') h1(`h1') lambda(`lambda') k(5) lag(4) vmat("nw")
```
Program can handle multiple instruments using standard Stata norms (B-splines are also used for these coefficients and the penalty matrix is appended accordingly). A graph will only be produced for the first endogenous variable listed unless the `mult` option is specified.  

## Inference 

Inference is not well-defined in light of cross-validation. SmoothLP by construction will give more narrow confidence intervals than LP, but as shown in recent work by [Montiel Olea, Qian,Plagborg-Møller, and Wolf](https://arxiv.org/abs/2503.17144) this is necessarily a consequence of inducing bias. More specifically, they show SmoothLP severely undercovers when following the heuristic reccomendations of BB19 (compute the variance-covariance matrix Newey-West w/ a penalization parameter that is 10% of what was selected during cross-validation). So an open question is how should standard errors be computed. 

To decide on default behavior, I ran simulations using the [GitHub repo](https://github.com/ckwolf92/lp_var_nberma/tree/main) of MOQPMW to see what choices could improve the coverage properties. To those ends, here are the default behaviors (and reccomendations)
  * The default is Huber-White standard errors. This should be paired with including a healthy number of lags as controls. 
  * The default undersmoothing is 1%; see below for a graph of how different undersmoothing choices affect coverage (what's plotted is the fraction of DGPs for with the SmoothLP estimator has proper coverage for different undersmoothing choices).  

These findings are consistent with broader results in the literature. Newey-West standard errors are the default in time series, but [Herbst and Johannsen (2024)](http://www.sciencedirect.com/science/article/pii/S0304407624000010) finds the NW variance matrix will often be biased while [Plagborg-Møller and Montiel Olea (2021)](https://joseluismontielolea.com/lp_inference_ecta.pdf) show that if a sufficient number of lags are included as controls, the usual HW errors are unbiased and have autocorrelation robustness properties. My simulations showed severe udnercoverage for NW relative to HW (for the same undersmoothing choices). 

## Syntax 

```
syntax anything(equalok) [if] [in], [Lambda(numlist) H(integer 20) K(integer 5) H1(integer 0) R(integer 2) Lag(integer 0) NWLag(integer 0) bdeg(integer 3) vmat(string) irfscale(integer 1) usmooth(real .01) ztail(real .05) MULT CUM NODRAW NOADJ]
```
* You can call `slp_irf` just like `reg`. To plot the IRF of `y` to `x`, list `y x` in that order. Every variable listed after `x` will be included in the list of controls. See section above for IV option
* H is the horizon length
* Lambda is a list of numbers considered for the ridge parameter $\lambda$ 
* K is the number of cross-validation folds
* H1 is the period the IRF starts
* r sets order of the limit polynomial. More specifically, the $r$-th derivitive of the IRF converges to 0 as $\lambda$ grows 
* Lag allows you to include lags of control variables in the conditioning set. For this, a `tsset` command must be run before `slp_irf`
* vmat takes option "nw" if you would rather use Newey-West standard errors over Huber-White, see the previous section for more discussion. 
  * You can specify Newey-West with `p` lags using `nwlag(p)`. Default is `H` as in the original R code.  
* irfscale allows you to scale the size of the shock (by default, it's a 1 std shock).
* adjstd calculates standard errors using a smaller penalty. Default is .1 (10% smaller) in line with BB19's reccomendations. This is because inference is not well-defined because of the cross-validation. 
* ztail allows for the confidence bands to be adjusted. Default is .05, which corresponds to 90% confidence intervals. 
* Alternativly, you can add the `noadj` option for a pure plot of the coefficients
* Add `cum` to instead do cumulative IRFs
* If you have multiple IVs and you want graphs for all endogenous variables, add `mult`
* Add `nodraw` for no graphs

To recap and list all of the things that are stored following the command 
* matrix results: the coeficients that are plotted
* matrix ic: the bands on the coeficients
* results and ic are also stored as variables, along with lagged variables if that option is chosen 
* matrix se: the standard errors
* scalar min_lambda: the $\lambda$ parameter chosen through cross-validation. Note the penalization is scaled by the number of observations across all horizon estimates (will be pretty close to $T * (H+1-h_1)$ ) in order to keep $\lambda$ relatively small. 

## Future Development 

 I will be updating it as I work through a couple projects. Some things I have planned 

* more customizability: confidence bands, graphs, ability to have different lag lengths across controls  
* error messages that are common practice in stata packages (e.g., if Lag() is specified but data not loaded with time series format) 

An important disclaimer for these programs is that it may be advisable to transform your variables if they are extremely small in magnitude due to well-known precision issues associated with Stata/Mata. This can happen in general (e.g., 7th digit of 1/6) but is particularly an issue when inverting matricies.  


***
Thank you to my RA, Claud, for all their help!

