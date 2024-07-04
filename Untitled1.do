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

timer on 1

lproj_size dlip `dynLag' (ggg3 ggg2 ggg4 ggg1 = gg3 gg2 gg4 gg1), h(24) lambda(0.1) k(5) mult noadj cum

timer off 1 

timer list 

*reg ggg1hat temp1 `dynLag' `extra' 

scalar ND = _se[$ffrSN]

In

putexcel set "C:\Users\pblit\OneDrive\Documents\huge_matrix.xlsx", sheet("M") replace
putexcel A1=matrix(X)


forvalues i = 0/`=hor' {

 foreach var of local outcome  {
 	forvalues j = 1/`=_N-`i'' {
    quietly sum `var' in `j'/`=`j'+`i''
    quietly replace temp = r(sum) in `j'
}
	display "`var'.`i'"
 	
     quietly ivregress gmm temp L(1/`pp').tffr `dynLag' (ggg3 ggg2 ggg4 ggg1 = gg3 gg2 gg4 gg1), vce(robust)
	 
	 matrix V = e(V)
	 
	 scalar ND = _se[$ffrSN]
	 scalar PD = _se[$ffrSP]

     scalar Nb`var'h`i' = (_b[$ffrBN]-_b[$ffrSN])/ND
     scalar Nse`var'h`i' = sqrt(V[b_n,b_n]+V[s_n,s_n] - 2 * V[s_n,b_n])/ND
  
	 scalar Pb`var'h`i' = (_b[$ffrBP]-_b[$ffrSP])/PD
	 scalar Pse`var'h`i' = sqrt(V[b_p,b_p]+V[s_p,s_p] - 2 * V[s_p,b_p])/PD
	 
	 quietly replace Nb`var' = Nb`var'h`i' if h==`i'
     quietly replace Nup90b`var' = Nb`var'+1.68*Nse`var'h`i' if h==`i'
     quietly replace Nlo90b`var' = Nb`var'- 1.68*Nse`var'h`i' if h==`i'

	 quietly replace Pb`var' = Pb`var'h`i' if h==`i'
	 quietly replace Pup90b`var' = Pb`var'+1.68*Pse`var'h`i' if h==`i'
	 quietly replace Plo90b`var' = Pb`var'- 1.68*Pse`var'h`i' if h==`i'
	 
	 scalar Nb`var'h`i' = (b * _b[$ffrBN]- a * _b[$ffrSN])
     scalar Nse`var'h`i' = sqrt(V[b_n,b_n] * b^2 + V[s_n,s_n] * a^2 - 2 * a * b * V[s_n,b_n])
  
	 scalar Pb`var'h`i' = (d * _b[$ffrBP]- c * _b[$ffrSP])
	 scalar Pse`var'h`i' = sqrt(V[b_p,b_p] * d^2 +V[s_p,s_p] * c^2 - 2 * c * d * V[s_p,b_p])
	 
	 quietly replace LNb`var' = Nb`var'h`i' if h==`i'
     quietly replace LNup90b`var' = LNb`var'+1.68*Nse`var'h`i' if h==`i'
     quietly replace LNlo90b`var' = LNb`var'- 1.68*Nse`var'h`i' if h==`i'

	 quietly replace LPb`var' = Pb`var'h`i' if h==`i'
	 quietly replace LPup90b`var' = LPb`var'+1.68*Pse`var'h`i' if h==`i'
	 quietly replace LPlo90b`var' = LPb`var'- 1.68*Pse`var'h`i' if h==`i'
  }
}

label variable dlip "Industrial Production"
label variable dlcpi "Inflation"
**label variable dgs1 "1Y Yield"
label variable debp "Excess Bond Premium"
label variable dcons "Consumption"

local outcome dlip dcons debp dlcpi

foreach var of local outcome {
    local varlabel : variable label `var' // Dynamic title using the variable name

    tw (rarea Nup90b`var' Nlo90b`var' h, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter Nb`var' h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if h<=hor, ///
         title(`varlabel') ///
         saving(Nvargk_`var'.gph, replace) nodraw 
		 
	tw (rarea Pup90b`var' Plo90b`var' h, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter Pb`var' h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if h<=hor, ///
         title(`varlabel') ///
         saving(Pvargk_`var'.gph, replace) nodraw
}


graph combine Nvargk_dcons.gph Nvargk_dlip.gph Nvargk_dlcpi.gph, cols(3) graphregion(fcolor(255 255 244)) title("Negative Shock Size Effect ({&sigma})") nodraw

// Save the combined "Negative" graph
graph save "Negative.gph", replace

// Combine the 4 graphs for the "Positive" set
graph combine Pvargk_dcons.gph Pvargk_dlip.gph Pvargk_dlcpi.gph, cols(3) graphregion(fcolor(255 255 244)) title("Positive Shock Size Effect ({&sigma})") nodraw

// Save the combined "Positive" graph
graph save "Positive.gph", replace

// Combine the two sets of graphs vertically
graph combine "Negative.gph" "Positive.gph", rows(2) cols(1) ysize(2.5) xsize(4) graphregion(fcolor(255 255 244))


foreach var of local outcome {
    local varlabel : variable label `var' // Dynamic title using the variable name

    tw (rarea LNup90b`var' LNlo90b`var' h, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter LNb`var' h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if h<=hor, ///
         title(`varlabel') ///
         saving(LN_`var'.gph, replace) nodraw 
		 
	tw (rarea LPup90b`var' LPlo90b`var' h, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter LPb`var' h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if h<=hor, ///
         title(`varlabel') ///
         saving(LP_`var'.gph, replace) nodraw
}


graph combine LN_dcons.gph LN_dlip.gph LN_dlcpi.gph, cols(3) graphregion(fcolor(255 255 244)) title("Negative Shock Size Effect (Level)") nodraw

// Save the combined "Negative" graph
graph save "LNegative.gph", replace

// Combine the 4 graphs for the "Positive" set
graph combine LP_dcons.gph LP_dlip.gph LP_dlcpi.gph, cols(3) graphregion(fcolor(255 255 244)) title("Positive Shock Size Effect (Level)") nodraw

// Save the combined "Positive" graph
graph save "LPositive.gph", replace

// Combine the two sets of graphs vertically
graph combine "LNegative.gph" "LPositive.gph", rows(2) cols(1) ysize(2.5) xsize(4) graphregion(fcolor(255 255 244))

capture log close
