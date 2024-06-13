program lproj_irf, eclass 
    version 15.0
    syntax anything(equalok) [if] [in], H(integer) Lambda(numlist) K(integer) [H1(integer 0)] [R(integer 2)] [Lag(integer 0) bdeg(integer 3) vmat(string) se(integer 1) MULT CUM NODRAW NOADJ]

    * Parse input
    local y : word 1 of `anything'
    // initialize x 
    local x 
    local H = `h'
    local h1 = `h1'
    // Remove the dependent variable from the variable list
    local mesh: subinstr local anything "`y'" "", word
    if (strpos("`mesh'", "(") > 0) {
    	local pos_open = strpos("`mesh'", "(")

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
    local lambda `lambda'
    local K = `k'
    local Lag = `lag'
    local w `contr'
    local vmat = "`vmat'"
    
    if (`Lag' > 0) {
    	foreach var of local varlist {
		forv i = 1/`Lag' {
		quietly gen `var'_`i' = L`i'.`var'
		local w `w' `var'_`i'	
		} 
	}
    }
    	
    
    scalar delt = 0 
    
    if (strpos("`mesh'", "(") > 0) {

 // isolate parentheses expression 
        local paren = substr("`mesh'",`pos_open', .) 
	local paren: subinstr local paren "(" ""
        local paren: subinstr local paren ")" ""
	// isolate endogenous and exogenous variables 
	local pos_eq = strpos("`paren'", "=")
	local endg = substr("`paren'",1, `=`pos_eq'-1') 
	local exg = substr("`paren'",`=`pos_eq'+1', .) 
	local EV = wordcount("`endg'")
	if ("`mult'" != "") {
		forvalues i=0/`=`EV' -1' {
			local yy : word `=`EV' -`i'' of `endg'
			local xx : word `=`EV' -`i'' of `exg' 
			quietly reg `yy' `xx' `w'
			scalar delt = e(rmse)
			quietly predict `yy'_hat
			local x `yy'_hat `x'
		}
	}
	else {
		local yy : word  1 of `endg'
		local xx : word 1 of `exg' 
		quietly reg `yy' `xx' `w'
		scalar delt = e(rmse)
		quietly predict `yy'_hat
		local x `yy'_hat 
		forvalues i=0/`=`EV' -2' {
			local yy : word `=`EV' -`i'' of `endg'
			local xx : word `=`EV' -`i'' of `exg' 
			quietly reg `yy' `xx' `w'
			quietly predict `yy'_hat
			local w `yy'_hat `w'
		}
		local EV  = 1 
	}
    }
    else {
    	local EV  = 1 
    }
    
    * Additional initializations
    
     //we need to drop all places we have missing values since we're doing matrix ops
    // so first we will create the y matrix to maximize data at longer horizons 
    
    quietly drop if _n <= `Lag'
    quietly drop if missing(`y')
    
    local T = _N
    
    if ("`cum'" != "") {
    	forvalues i=`h1'/`H' {
		quietly gen temp_`i' = .
		forvalues j=1/`=`T'-`i'' {
			quietly sum `y' in `j'/`=`j'+`i''
			quietly replace temp_`i' = r(sum) in `j'
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
    
    local bdeg = `bdeg'

    * Create the range of values for h
    tempvar h_range
    quietly generate `h_range' = `h1' + _n - 1 if _n <= `H'+1-`h1'

    * Generate basis functions
    quietly bspline, xvar(`h_range') power(`bdeg') knots(`=`h1''(1)`=`H'') gen(bs)
    
    mkmat bs*, matrix(bs)
    

    * Create a matrix from the basis variables
    matrix basis = bs[1..`=`H'+1-`h1'',`=2-`h1''...]
    
    drop bs*
    * Create covariates matrix if w is specified
    matrix w = J(_N, 1, 1) 
    
        if ("`w'" != "") {
        tempname wmat
        mkmat `w', matrix(`wmat')
        matrix w = (J(_N, 1, 1) , `wmat')
    }
    else {
        matrix w = J(_N, 1, 1)
    }
    
    * 1 shock std dev
    if (("`w'" != "") & (strpos("`mesh'", "(")==0)) {
        quietly reg `x' `w'
        scalar delt = e(rmse)
    }
    else {
        if ("`w'" == "") {
		scalar delt = r(sd)
	}
    }
    
    if ("`noadj'" != ""){
    	local delta = 1
    }
    else {
    	local delta = `=delt'*`se'
    }
    	

	* Additional initializations
    local T = _N
    local HR = `H' + 1 - `h1'
    local TS = `T' * `HR'
    local XS = colsof(basis)*`EV'
    local NW = colsof(w)
    
    local back = `T'-`h1'
    local L = wordcount("`lambda'")
    mkmat `x', matrix(x)
    
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
    mata: cvtwirl(`T',Y,X,P,basis,IDX,`h1',`H',`L',`K',`TS',`XS',`delta',`EV',lambda_vec, "`vmat'","`mult'")
    
        * Prepare data for graphing
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
         (scatter result1 `sc_time', c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if (`shift'<=time)&(time<=`=`i'*(`H'+1)-1'), title("IRF of `y' for shock to `xx'") name("`xx'") xtitle(horizon)
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
	IDX = J(TS, 2, 0)
	Y = J(TS, 1, 0)
	Xb = J(TS, XS, 0)
	Xc = J(HR,width,0)
	II = I(HR)
	
		
	for(t=1; t<=back; t++) {
		idx_beg = (t-1)*HR + 1
		idx_end = t*HR
    
		IDX[|idx_beg,1 \ idx_end,1|] = J(HR, 1, t)
		IDX[|idx_beg,2 \ idx_end,2|] = range(h1, H,1)
   
		Y[|idx_beg,1 \ idx_end,1|] = yy[|t,1 \ t,HR|]'
		
		stack = basis*x[t,1]
		for(i=2; i<=EV; i++) {
			stack = stack, basis*x[t,i]
		}
    
		Xb[|idx_beg,1 \ idx_end,XS|] = stack
		for(i=1; i<=NW; i++){
			cdx_beg = (i-1)*HR + 1
			cdx_end = i*HR
			Xc[|idx_beg,cdx_beg \ idx_end,cdx_end|] = II * w[t,i]
		}
		Xxc = J(HR,width,0)
		Xc = Xc\Xxc
	}
	X = Xb, Xc[|1,1 \ TS,width|] 
	st_matrix("X", X)
	st_matrix("Y", Y) 
	st_matrix("IDX", IDX)
	XX = cross(X, X)
	XY = cross(X, Y) 
	P = J(cols(X), cols(X), 0)
        D = I(XS)
        for (i=1; i<=r; i++) {
		D = D[2..rows(D), .] - D[1..rows(D)-1, .]
	}
	DD = D' * D
        P[1::XS, 1::XS] = DD
	st_matrix("P",P)
} 



void function cvtwirl( real scalar T,
                     real matrix Y,
		     real matrix X,
		     real matrix P,
		     real matrix basis,
		     real matrix IDX,
		     real scalar h1,
		     real scalar H, 
		     real scalar L,
		     real scalar K,
		     real scalar TS,
		     real scalar XS,
		     real scalar delta,
		     real scalar EV,
		     string vector lambda_vec,
		     string vmat,
		     string mult)
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
			A = cross(X_in, X_in) + lambda_l * rows(Y_in) * ((K - 1) / K) * P
			b = cross(X_in, Y_in)
			beta = invsym(A) * b
			rss_l[i,1] = mean((Y_out - X_out * beta):^2)
		}
		avg_rss = mean(rss_l)
		//avg_rss
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
        lambda_opt = min_lambda 

	XX = cross(X, X)
	XY = cross(X, Y)
        A = XX + lambda_opt * rows(Y) * P
        b = XY
        theta = invsym(A) * b
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
	

        u = Y - X * theta
        S = X :* (u * J(1, cols(X), 1))

        bread = invsym(XX+ lambda_opt * rows(Y) * P)

        nlag = H
	lagseq = 0::nlag
        V = cross(S, S)
	meat = V 
	if (vmat=="nw"){
		lagseq = 0::nlag
		weights = 1 :- lagseq :/ (nlag+1)
		for (i=1; i<=nlag; i++){
		Gammai = cross(S[(i+1)::rows(S),], S[1::(rows(S)-i),])
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
        conf[,1] = mu :+ se * invnormal(0.10)
        conf[,2] = mu :+ se * invnormal(0.90)
	

        irc = J(rows(se)+1, 2, .)
        irc[(h1+1)::linked,] = conf * delta
	results[|(h1+1),1 \ linked,1 |] = mu * delta

        st_matrix("results", results)
        st_matrix("irc", irc)
	st_matrix("se", se)
        
}

end 
