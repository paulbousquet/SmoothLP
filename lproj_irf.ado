program lproj_irf, eclass 
    version 15.0
    syntax anything(equalok) [if] [in], H(integer) Lambda(numlist) K(integer) [H1(integer 0)] [R(integer 2)] [Lag(integer 0) bdeg(integer 3) vmat(string)]

    * Parse input
    local y : word 1 of `anything'
    local H = `h'
    local h1 = `h1'
    // Remove the dependent variable from the variable list
    local mesh: subinstr local anything "`y'" "", word
    if (strpos("`mesh'", "(") > 0) {
    	local pos_open = strpos("`mesh'", "(")

 // Extract the control variables
        local contr = substr("`mesh'", 1, `pos_open' - 1)

	local contr `contr' `y'
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
	local yy : word  1 of `endg'
	local xx : word 1 of `exg' 
	quietly reg `yy' `xx' `w'
	scalar delt = e(rmse)
	predict `yy'_hat
	local x `yy'_hat 
	
	forvalues i=0/`=`EV' -2' {
		local yy : word `=`EV' -`i'' of `endg'
		local xx : word `=`EV' -`i'' of `exg' 
		quietly reg `yy' `xx' `w'
		predict `yy'_hat
		local w `yy'_hat `w'
	}
    }
    else {
    	local EV  = 1 
    }
    
    local shabang `x' `w' `y'
    
    quietly drop if _n <= `Lag'
     
    foreach var of local shabang {
    	quietly drop if missing(`var')
    }
    
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
    
    local delta = `=delt'

    * Additional initializations
    local HR = `H' + 1 - `h1'
    local TS = _N * `HR'
    local XS = colsof(basis)
    local NW = colsof(w)
    local T = _N
    local back = `T'-`h1'
    local L = wordcount("`lambda'")
    mkmat `y', matrix(y)
    mkmat `x', matrix(x)
    
    mata: y = st_matrix("y")
    mata: x = st_matrix("x")
    mata: x = x[1...,1]
    mata: w = st_matrix("w")
    mata: basis = st_matrix("basis")
    mata: lambda_vec = tokens("`lambda'")
    mata: twirl(`back',y,x,w,basis,`h1',`H',`HR',`NW',`TS',`XS',`r')
    mata: X = st_matrix("X")
    mata: Y = st_matrix("Y") 
    mata: P = st_matrix("P") 
    mata: IDX = st_matrix("IDX")
    mata: cvtwirl(`T',Y,X,P,basis,IDX,`h1',`H',`L',`K',`TS',`XS',`delta',lambda_vec, "`vmat'")
    
        * Prepare data for graphing
    svmat double results, names(result)
    svmat double irc, names(irc)
    gen time = _n - 1


       
         tw (rarea irc1 irc2 time, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
         (scatter result1 time, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick) legend(off) graphregion(fcolor(255 255 244))) if time<=`H'

end



mata: 

void function twirl( real scalar back,
                     real matrix y,
		     real matrix x,
		     real matrix w,
		     real matrix basis,
		     real scalar h1,
		     real scalar H, 
		     real scalar HR,
		     real scalar NW,
		     real scalar TS,
		     real scalar XS,
		     real scalar r)
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
    
   
		y_range_start = t+h1
		y_range_end = min((t+H, rows(y)))
		Y[|idx_beg,1 \ idx_end,1|] = (y[|y_range_start,1 \ y_range_end,1|] \ J(HR-(			y_range_end-y_range_start+1), 1, 0))
    
		Xb[|idx_beg,1 \ idx_end,XS|] = basis*x[t]
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
		     string vector lambda_vec,
		     string vmat)
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
	linked = H + 1
        results = J(linked, 1, 0)
        theta = J(cols(X), 1, 0)
	ml = st_numscalar("min_lambda")
        lambda_opt = ml

	XX = cross(X, X)
	XY = cross(X, Y)
        A = XX + lambda_opt * rows(Y) * P
        b = XY
        theta = invsym(A) * b
        beta = theta[1..XS, 1]
	mu = st_matrix("basis") * beta 
        

        u = Y - X * theta
        S = X :* (u * J(1, cols(X), 1))

        bread = invsym(XX+ lambda_opt * rows(Y) * P)

        nlag = H
	lagseq = 0::nlag
        V = cross(S, S)
	meat = V 
	fok = 5
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

        V = basis * V[1::XS, 1::XS] * basis'
        se = sqrt(diagonal(V))

        conf = J(rows(se), 2, .)
        conf[,1] = mu :+ se * invnormal(0.10)
        conf[,2] = mu :+ se * invnormal(0.90)

        irc = J(H+1, 2, .)
        irc[(h1+1)::(H+1),] = conf * delta
	
	results[|(h1+1),1 \ linked,1 |] = mu * delta

        st_matrix("results", results)
        st_matrix("irc", irc)
	st_matrix("se", se)
	st_matrix("theta",theta)
        
}

end 
