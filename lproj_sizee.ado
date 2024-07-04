program lproj_sizee, eclass 
    version 15.0
    syntax anything(equalok) [if] [in], [H(integer 24) K(integer 5) H1(integer 0) R(integer 2) bdeg(integer 3)]

    * Parse input
	
    local y : word 1 of `anything'
    // initialize x 
    local x 
    local trc 
    local H = `h'
    local h1 = `h1'
    // Remove the dependent variable from the variable list
    local mesh: subinstr local anything "`y'" "", word
	local ivdum = strpos("`mesh'", "(") 
    if (`ivdum'> 0) {
    	local pos_open = `ivdum'

 // Extract the control variables
        local contr = substr("`mesh'", 1, `pos_open' - 1)

	local contr `contr'
    }
    else {
    	local x : word 1 of `mesh'
	local contr : subinstr local mesh "`x'" "", word
    }

    local varlist `contr' `x' `y'
    
    if (`h1'> 0) {
    	local contr `contr' `y'
    }
    
    
    
    local r = `r'
    local K = `k'
    local w `contr'
    
    	
    
    scalar delt = 0 
	


 // isolate parentheses expression 
        local paren = substr("`mesh'",`pos_open', .) 
	local paren: subinstr local paren "(" ""
        local paren: subinstr local paren ")" ""
	// isolate endogenous and exogenous variables 
	local pos_eq = strpos("`paren'", "=")
	local endg = substr("`paren'",1, `=`pos_eq'-1') 
	local exg = substr("`paren'",`=`pos_eq'+1', .) 
	local EV = wordcount("`endg'")
		forvalues i=1/`EV' {
			local yy : word `i' of `endg'
			quietly reg `yy' `exg' `w'
			scalar delt = e(rmse)
			quietly predict hat_`yy'
			local x `x' hat_`yy'
			local trc `trc' `yy'
		} 
    
    * Additional initializations
    
     //we need to drop all places we have missing values since we're doing matrix ops
    // so first we will create the y matrix to maximize data at longer horizons 
    
    quietly drop if missing(`y')
    local T = _N

    	forvalues i=`h1'/`H' {
		quietly gen temp_`i' = .
		forvalues j=1/`=`T'-`i'' {
			quietly sum `y' in `j'/`=`j'+`i''
			quietly replace temp_`i' = r(sum) in `j'
		}
		}
    
    local shabang `x' `w' 
     
    foreach var of local shabang {
    	quietly drop if missing(`var')
    }
	
    marksample touse 
    
    tempvar esample 
    gen byte `esample' = `touse'

    
    mkmat temp_*, matrix(yy)
	drop temp_*

    * Create the range of values for h
    tempvar h_range
    quietly generate `h_range' = `h1' + _n - 1 if _n <= `H'+1-`h1'

    * Generate basis functions
    quietly bspline, xvar(`h_range') power(`bdeg') knots(`=`h1''(1)`=`H'') gen(bs)
    
    mkmat bs*, matrix(bs)
    

    * Create a matrix from the basis variables
    matrix basis = bs[1..`=`H'+1-`h1'',1'...]
	    
    drop bs*
    * Create covariates matrix if w is specified
    matrix w = J(_N, 1, 1) 
    
        tempname wmat
        mkmat `w', matrix(`wmat')
        matrix w = (J(_N, 1, 1) , `wmat')
   
    
	local delta = 1 
	
	 
    	

	* Additional initializations
    local T = _N
    local HR = `H' + 1 - `h1'
    local TS = `T' * `HR'
    local XS = colsof(basis)*`EV'
    scalar XSt = colsof(basis)
    local NW = colsof(w)
    
    local back = `T'-`h1'
    mkmat `x', matrix(x)

    mkmat `trc', matrix(xz)

	
    
    mata: yy = st_matrix("yy")
    mata: x = st_matrix("x")
    mata: w = st_matrix("w")
    mata: basis = st_matrix("basis")
    mata: twirl(`back',yy,x,w,basis,`h1',`H',`HR',`NW',`TS',`XS',`r',`EV')
    tempname bb
    ereturn post, esample(`esample')
    matrix `bb' = beta'
    ereturn post `bb'
    ereturn local cmd="bootstrap"
    drop hat_*
end



mata: 

real scalar function idb(idx,blk){
	return ((idx-1)*blk+1)
}


void function twirl( real scalar back,
                     real matrix yy,
		     real matrix x,
		     real matrix w,
		     real matrix basis,
		     real scalar h1,
		     real scalar H, 
		     real scalar HR,
		     real scalar NW,
		     real scalar TS,
		     real scalar XS,
		     real scalar r,
		     real scalar EV)
{
	width = NW * HR	
	IDX = J(TS, 2, .)
	Y = J(TS, 1, .)
	Xb = J(TS, XS, .)
	Xc = J(TS,width,.)
	II = I(HR)
	XSt = XS/EV
	
		
	for(t=1; t<=back; t++) {
		stt = (t-1)*HR + 1
		edd = t*HR
    
		IDX[|stt,1 \ edd,2|] = J(HR, 1, t), range(h1, H,1)
   
		Y[|stt,1 \ edd,1|] = yy[|t,1 \ t,HR|]'
		
		for(i=1; i<=EV; i++) {
			Xb[|stt,idb(i,XSt) \ edd,i*XSt|] = basis*x[t,i]
		}
	
		for(i=1; i<=NW; i++){
			Xc[|stt,idb(i,HR) \ edd,i*HR|] = II * w[t,i]
		}
	}

	X = Xb, Xc
	sel = Y :!= .
	TS = length(Y)
	P = J(cols(X), cols(X), 0)
        D = I(XSt)
        for (i=1; i<=r; i++) {
		D = D[2..rows(D), .] - D[1..rows(D)-1, .]
	}
	DD = D' * D
    P[1::XSt, 1::XSt] = DD
	for(i=2; i<=EV; i++) {
		stt = XSt*(i-1)+1
		edd = XSt*i 
		P[stt::edd,stt::edd] = DD 
	}
	XX = cross(X, X)
	XY = cross(X, Y)
	A = luinv(XX + .1 * rows(Y) * P)
	theta = A * XY
	beta = basis*theta[1..XSt,1]
	st_matrix("beta", beta)
	beta 
	
} 





end 
