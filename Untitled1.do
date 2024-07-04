drop _all
clear all


/*******************************************************************************
** RAW DATA IMPORTATION AND DATA SETUP

** SEE README SHEET OF EXCEL FILE FOR VARIABLE DEFINITIONS
*******************************************************************************/

import excel Monetarydat.xlsx, sheet("Monthly") firstrow case(l)

gen mdate = m(1959m1) + _n-1
tsset mdate, m

global shockN nrr
global shockP rr


gen ind1 = (tffrch>.1) 

gen ind2 = (tffrch<-.1)

gen ind3 = (-.1 <= tffrch & tffrch<= - .01) 

gen ind4 = (.1 >= tffrch & tffrch >= .01)

replace tffrch = tffrch * 100

gen gg1 = ind1 * $shockP

gen gg2 = ind2 * $shockN 

gen gg3 = ind3 * nswan

gen gg4 = ind4 * rr

gen gg5 = (1-ind3) * (1-ind4) * $shockP 

global shockBN gg2 

global shockSN gg3

global shockBP gg1

global shockSP gg4

global Oshock gg5

gen ggg1 = ind1 * tffr

gen ggg2 = ind2 * tffr

gen ggg3 = ind3 * tffr 

gen ggg4 = ind4 * tffr

global ffrBN ggg2 

global ffrSN ggg3

global ffrBP ggg1

global ffrSP ggg4

gen zlb = (tffr < .2)

gen t = _n

/*******************************************************************************
**  SET UP FOR LATER IRF PLOTTING
*******************************************************************************/
**local biascor 

local set lip lcpi cons ebp

local outcome 

foreach var of local set {
	gen d`var' = D.`var'
	local outcome `outcome' d`var'
}


local p = 1
local q = 1 /* q = 0 corresponds to recursiveness assumption, q = 1 is less restrictive */
local lg = 11
local pp = `lg' + `p'

local shocks $shockSN $shockBN $shockSP $shockBP

local ffrs $ffrSN $ffrBN $ffrSP $ffrBP
 
**gen dlip = 100 * D.lip
 
**local varcontemp : subinstr local outcome "lip" "", all
**local varcontemp `varcontemp' dlip

local extra unemp zlb sent mpu

local varcontemp `extra' 

local varlag `outcome' tffr

local allvar `set' `extra'

* Initialize an empty list for storing dynamic lag terms
local dynLag

* Iterate over each var in allvar to create dynamic lag terms
foreach var of local varlag {
	forvalues i=`q'/`pp' {
		quietly gen `var'_`i' = L`i'.`var'
		local dynLag `dynLag' `var'_`i'
	}
}

foreach var of local varcontemp {
	forvalues i=0/`pp' {
		quietly gen `var'_`i' = L`i'.`var'
		local dynLag `dynLag' `var'_`i'
	}
}
	

foreach var of local set {
	forvalues i=`q'/`p' {
		quietly gen `var'_`i' = L`i'.`var'
		local dynLag `dynLag' `var'_`i'
	}
}

reg ggg1 gg1 `dynLag'

predict ggg1hat

drop if mdate < m(1973m3)


foreach var of local outcome  {

  quietly gen Nb`var' = .
  quietly gen Nup90b`var' = .
  quietly gen Nlo90b`var' = .

  
  quietly gen Pb`var' = .
  quietly gen Pup90b`var' = .
  quietly gen Plo90b`var' = .

  
  quietly gen LNb`var' = .
  quietly gen LNup90b`var' = .
  quietly gen LNlo90b`var' = .
  
  quietly gen LPb`var' = .
  quietly gen LPup90b`var' = .
  quietly gen LPlo90b`var' = .
}

scalar s_n = 1
scalar b_n = 2
scalar s_p = 3
scalar b_p = 4

scalar hor = 24

generate h = min(_n, hor)-1 

gen temp = . 

 	forvalues j = 1/`=_N' {
    quietly sum ffr in `j'/`=`j''
    quietly replace temp = r(sum) in `j'
}

ivregress gmm temp `dynLag' (ggg3 ggg2 ggg4 ggg1 = gg3 gg2 gg4 gg1), vce(robust)

scalar a = -.25 / _b[$ffrSN]

scalar b = -.5 / _b[$ffrBN]

scalar c = .25 / _b[$ffrSP]

scalar d = .5 / _b[$ffrBP]

local lambda 0.00001 0.0001 0.001 0.003 0.005 0.007 0.008 .009 .01 .02 .03 .06 

local lambda 0

set trace off

quietly gen temp1 = dlip*100

set trace off


*lproj_sizee dlip `dynLag' (ggg3 ggg2 ggg4 ggg1 = gg3 gg2 gg4 gg1)


* Step 4: Bootstrap the weighted average adjusted difference
bootstrap _b, noisily reps(50): lproj_sizee dlip `dynLag' (ggg3 ggg2 ggg4 ggg1 = gg3 gg2 gg4 gg1)
