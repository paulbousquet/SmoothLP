program slp_irf, eclass 
    version 15.0
    syntax anything(equalok) [if] [in], [Lambda(numlist) H(integer 20) K(integer 5) H1(integer 0) R(integer 2) Lag(integer 0) NWLag(integer 0) bdeg(integer 3) vmat(string) irfscale(integer 1) usmooth(real .01) ztail(real .05) MULT CUM NODRAW NOADJ]

    // Check if data is time series
    capture tsset
    if (_rc != 0 & `lag' > 0) {
        di as error "data must be tsset with lag option"
        exit 459
    }

    // Input validation
    if `h' <= 0 {
        di as error "h() must be a positive integer"
        exit 198
    }
    if `h1' < 0 {
        di as error "h1() must be non-negative"
        exit 198
    }
    if `h1' > `h' {
        di as error "h1() cannot be greater than h()"
        exit 198
    }
    if `ztail' <= 0 | `ztail' >= 1 {
        di as error "ztail() must be between 0 and 1"
        exit 198
    }
    if `k' <= 1 {
        di as error "k() must be greater than 1 for cross-validation"
        exit 198
    }
	if `usmooth' < 0 {
        di as error "adjustment must be non-negative"
        exit 198
    }

    // Parse input
    local y : word 1 of `anything'
    
    // Confirm dependent variable exists
    capture confirm variable `y'
    if _rc != 0 {
        di as error "variable `y' not found"
        exit 111
    }
    
    // Initialize x 
    local x 
	local trc 
    local H = `h'
    local h1 = `h1'
	local nlag = `h'
	
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
    
    // Confirm all variables exist
    foreach var of local varlist {
        capture confirm variable `var'
        if _rc != 0 {
            di as error "variable `var' not found"
            exit 111
        }
    }
    
    if (`h1'> 0) {
    	local contr `contr' `y'
    }
	
	if (`nwlag'>0){
		local nlag = `nwlag'
		
	}
	
    local r = `r'
    local lambda `lambda'
    local K = `k'
    local Lag = `lag'
    local w `contr'
    local vmat = "`vmat'"

    if "`lambda'" == "" {
		local lambda = 0
	}
    
    if (`Lag' > 0) {
    	foreach var of local varlist {
		forv i = 1/`Lag' {
		quietly gen `var'_`i' = L`i'.`var'
		local w `w' `var'_`i'	
		} 
	}
    }
    	
    scalar delt = 0 
	
    if (`ivdum' > 0) {

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
			quietly predict `yy'_hat
			local x `x' `yy'_hat 
			local trc `trc' `yy'
		}
    }
    else {
    	local EV  = 1 
    }
		 
    // Check for sufficient observations
    quietly count if !missing(`y')
    if r(N) <= `h' {
        di as error "insufficient observations: need at least `=`h'+1' non-missing observations"
        exit 2001
    }
    
     //we need to drop all places we have missing values since we're doing matrix ops
    // so first we will create the y matrix to maximize data at longer horizons 
    
    quietly drop if _n <= `Lag'
    local T = _N
    
    if ("`cum'" != "") {
    	forvalues i=`h1'/`H' {
		quietly gen temp_`i' = .
		forvalues j=1/`=`T'-`i'' {
			quietly sum `y' in `j'/`=`j'+`i''
			quietly replace temp_`i' = r(sum) / (r(N) == `i'+1) in `j'
		}
		}
	}
    else {
    	forvalues i=`h1'/`H' {
		quietly gen temp_`i'= F`i'.`y'
	}
    	
    }
    
    local shabang `x' `w' 
     
    foreach var of local shabang {
    	quietly drop if missing(`var')
    }
	
    mkmat temp_*, matrix(yy)
	drop temp_*

    // Create the range of values for h
    tempvar h_range
    quietly generate `h_range' = `h1' + _n - 1 if _n <= `H'+1-`h1'

    // Generate basis functions
    quietly bspline, xvar(`h_range') power(`bdeg') knots(`=`h1''(1)`=`H'') gen(bs)
    
    mkmat bs*, matrix(bs)

    // Create a matrix from the basis variables
    matrix basis = bs[1..`=`H'+1-`h1'',1'...]
	    
    drop bs*
    
    // Create covariates matrix if w is specified
    matrix w = J(_N, 1, 1) 
    
        if ("`w'" != "") {
        tempname wmat
        mkmat `w', matrix(`wmat')
        matrix w = (J(_N, 1, 1) , `wmat')
    }
    else {
        matrix w = J(_N, 1, 1)
    }
    
    // 1 shock std dev
    if (("`w'" != "") & (`ivdum'==0)) {
        quietly reg `x' `w'
        scalar delt = e(rmse)
    }
    else {
        if ("`w'" == "") {
		qui sum `x'
		scalar delt = r(sd)
	}
    }
    
    if ("`noadj'" != ""){
    	local delta = 1
    }
    else {
    	local delta = `=delt'*`irfscale'
    }
	
    // Additional initializations
    local T = _N
    local HR = `H' + 1 - `h1'
    local TS = `T' * `HR'
    local XS = colsof(basis)*`EV'
    local NW = colsof(w)
    
    local back = `T'-`h1'
    local L = wordcount("`lambda'")
    mkmat `x', matrix(x)
	if (`ivdum' > 0){
		mkmat `trc', matrix(xz)
	}
	else {
		matrix xz = (1,1)
	}
	
    mata: yy = st_matrix("yy")
    mata: x = st_matrix("x")
    mata: w = st_matrix("w")
    mata: basis = st_matrix("basis")
    mata: lambda_vec = tokens("`lambda'")
    mata: twirl(`back',yy,x,w,basis,`h1',`H',`HR',`NW',`TS',`XS',`r',`EV')
	mata: X = st_matrix("X")
    mata: Y = st_matrix("Y") 
    mata: P = st_matrix("P") 
    mata: IDX = st_matrix("IDX")
	mata: xz = st_matrix("xz")
	mata: sel = st_matrix("sel")
	mata: ivtwirl(`ivdum',`back',xz,basis,X,sel,`TS',`XS',`HR',`EV')
	mata: ZX = st_matrix("ZX")
	di "Data processing complete"
    mata: cvtwirl(`T',Y,X,P,basis,IDX,ZX,`h1',`H',`L',`K',`=TS',`XS',`delta',`EV',`nlag',`ztail',lambda_vec, "`vmat'","`mult'", `usmooth')
    
        // Prepare data for graphing
    svmat double results, names(result)
    svmat double irc, names(irc)
    gen time = _n - 1
    
    if ("`nodraw'" == ""){	
    	if ("`mult'" == ""){
		tw (rarea irc1 irc2 time, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter result1 time, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if time<=`H', title("IRF of `y' for shock to `x'") xtitle("horizon")
	}
	else{
		forvalues i=1/`EV' {
			local xx : word `i' of `endg'
			local shift = (`H'+1)*(`i'-1)
			tempvar sc_time 
			quietly gen `sc_time' = time - `shift'
		tw (rarea irc1 irc2 `sc_time', fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter result1 `sc_time', c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 255))) if (`shift'<=time)&(time<=`=`i'*(`H'+1)-1'), title("IRF of `xx' for shock to `xx'") name("`xx'") xtitle(horizon)
		}
	}
    }
        
end

mata: 

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
	IDX = select(IDX, sel)
	Y = select(Y, sel)
	X = select(X, sel)
	TS = length(Y)
	st_matrix("X", X)
	st_matrix("Y", Y) 
	st_matrix("IDX", IDX)
	st_numscalar("TS",TS)
	st_matrix("sel",sel)
	P = J(cols(X), cols(X), 0)
        D = I(XSt)
        for (i=1; i<=r; i++) {
		D = D[2..rows(D), 1..cols(D)] - D[1..rows(D)-1, 1..cols(D)]
	}
	DD = D' * D
    P[1::XSt, 1::XSt] = DD
	for(i=2; i<=EV; i++) {
		stt = XSt*(i-1)+1
		edd = XSt*i 
		P[stt::edd,stt::edd] = DD 
	}
	st_matrix("P",P)
} 

void function cvtwirl( real scalar T,
                     real matrix Y,
		     real matrix X,
		     real matrix P,
		     real matrix basis,
		     real matrix IDX,
			 real matrix ZX, 
		     real scalar h1,
		     real scalar H, 
		     real scalar L,
		     real scalar K,
		     real scalar TS,
		     real scalar XS,
		     real scalar delta,
		     real scalar EV,
			 real scalar nlag,
			 real scalar ztail,
		     string vector lambda_vec,
		     string vmat,
		     string mult,
			 real scalar usmooth)
{
	min_rss = .
	min_lambda = .
	ind = ceil((IDX[|1,1 \ TS,1|] :/ T) :* K)
	for (l=1; l <= L; l++){
		lambda_l = strtoreal(lambda_vec[l])
		rss_l = J(K, 1, .)
		for (i=1; i <= K; i++){
			mask_in = (ind :!= i)
			mask_out = (ind :== i)
			Y_in = select(Y, mask_in)
			X_in = select(X, mask_in)
			Y_out = select(Y, mask_out)
			X_out = select(X, mask_out)
			A = quadcross(X_in, X_in) + lambda_l * rows(Y_in) * ((K - 1) / K) * P
			b = quadcross(X_in, Y_in)
			beta = lusolve(A,b)
			rss_l[i,1] = mean((Y_out - X_out * beta):^2)
		}
		avg_rss = mean(rss_l)
		if ((min_rss == .) | (avg_rss < min_rss)) {
			min_rss = avg_rss
			min_lambda = lambda_l
		}
	}
	st_numscalar("min_rss", min_rss)
	st_numscalar("min_lambda", min_lambda)
	linked = (H + 1)*EV
        results = J(linked, 1, 0)
        theta = J(cols(X), 1, 0)
        lambda_opt = usmooth* min_lambda + 10^(-10)
	
	XX = quadcross(X, X)
	XY = quadcross(X, Y)
	
	
	
		
        A = XX + lambda_opt * rows(Y) * P
		bread = luinv(A)
        theta = bread * XY
		thetav = theta 
	beta = theta[1..XS, 1]
	if (mult != ""){
		mu = basis * beta[1..XS/EV,1] 
		for (i=2; i<=EV; i++){
			bigbas = basis * beta[((i-1)*XS/EV+1)..(i*XS/EV),1]
			mu = mu \ bigbas
		}
	}
	else{
		mu = basis * beta
	}
	

        u = Y - ZX * thetav
        S = X :* (u * J(1, cols(X), 1))

	lagseq = 0::nlag
        V = quadcross(S, S)
	meat = V 
	if (vmat=="nw"){
		lagseq = 0::nlag
		weights = 1 :- lagseq :/ (nlag+1)
		for (i=1; i<=nlag; i++){
		Gammai = quadcross(S[(i+1)::rows(S),], S[1::(rows(S)-i),])
		GplusGprime = Gammai + Gammai'	
		V = V + weights[i+1] * GplusGprime
		}
		meat = V
	}
	
	
        V = bread * meat * bread
	
	
	if (mult != ""){
		see = sqrt(diagonal(basis * V[1..XS/EV,1..XS/EV] * basis'))
		for (i=2; i<=EV; i++){
			start = (i-1)*XS/EV+1
			endd = i*XS/EV
			bigbas = sqrt(diagonal(basis * V[start..endd,start..endd] * basis')) 
			see = see \ bigbas
		}
		se = see
	}
	else{
		V = basis * V[1::XS, 1::XS] * basis'
		 se = sqrt(diagonal(V))
	}
	

       

        conf = J(rows(se), 2, .)
        conf[,1] = mu :+ se * invnormal(ztail)
        conf[,2] = mu :+ se * invnormal(1-ztail)
	

        irc = J(rows(se)+1, 2, .)
        irc[(h1+1)::linked,] = conf * delta
	results[|(h1+1),1 \ linked,1 |] = mu * delta
		
        st_matrix("results", results)
        st_matrix("irc", irc)
	st_matrix("se", se)
        
}

void function ivtwirl( real scalar ivdum,
                     real scalar back,
		     real matrix xz,
		     real matrix basis,
			 real matrix X, 
			 real matrix sel,
			 real scalar TS, 
		     real scalar XS,
			 real scalar HR, 
		     real scalar EV) 
{
	if (ivdum > 0){
		Xb = J(TS, XS, 0)
		for(t=1; t<= back; t++){
			idx_beg = (t-1)*HR + 1
			idx_end = t*HR
			stack = basis*xz[t,1]
			for(i=2; i<=EV; i++){
				stack = stack, basis*xz[t,i]
			}
			Xb[|idx_beg,1 \ idx_end,XS|] = stack
		}
		Xb = select(Xb,sel)
		ZX = Xb, X[1..rows(X),(XS+1)..cols(X)]
		st_matrix("ZX", ZX)
	}
	else{
	   ZX = X	
	}
	st_matrix("ZX",ZX)
}

real scalar function idb(idx,blk){
	return ((idx-1)*blk+1)
}
end 
