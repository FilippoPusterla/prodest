*! version 1.0.1 27Sep2016
*! version 1.0.2 05Jun2017 major changes in the code and how the whole routine works, added exponential and parameters options
*! authors: Gabriele Rovigatti, University of Chicago Booth, Chicago, IL & EIEF, Rome, Italy. mailto: gabriele.rovigatti@gmail.com
*!          Vincenzo Mollisi, Bolzano University, Bolzano, Italy & Tor Vergata University, Rome, Italy. mailto: vincenzo.mollisi@gmail.com

/***************************************************************************
** Stata program for prodest Postestimation
**
** Programmed by: Gabriele Rovigatti
**************************************************************************/
cap program drop prodest_p
program define prodest_p, sortpreserve eclass

	version 10.0

	syntax [anything] [if] [in] [, 			   ///
		RESIDuals 				   /// 
		EXPonential 			   ///
		PARameters					///
		]

	marksample touse 	// this is not e(sample)
	
	tempvar esample
	qui gen byte `esample' = e(sample)
	
	loc varlist `anything'
	
	loc mod = "`e(PFtype)'"
	
	/* check for options in launching command */
	if ( "`residuals'" == ""  & "`exponential'" == "" & "`parameters'" == ""){
		di as error "You must specify RESIDuals, EXPonential or PARameters"
		exit 198
	}
	/* check for correct usage of options */
	if ( ("`residuals'" != "" | "`exponential'" != "") & "`parameters'" != ""){
		di as error "the 'parameters' option cannot be used with other options"
		exit 198
	}
	
	if "`mod'" == "Cobb-Douglas"{ /* PART I: COBB-DOUGLAS */
		if ("`residuals'" != "" | "`exponential'" != "") {
			tempname beta
			mat `beta' = e(b)
			tempvar rhs 
			mat score double `rhs' = `beta'
			loc lhs `e(depvar)'
			if "`exponential'" != ""{
				qui gen `varlist' = exp(`lhs' - `rhs') `if'
			}
			else{
				qui gen `varlist' = `lhs' - `rhs' `if'
			}
		}
		else{ /* 'parameters' with cobb-douglas PF yields the results' table */
			di _coef_table, level($S_level)
		}
	}
	else { /* PART II: TRANSLOG */
		loc free = "`e(free)'"
		loc state = "`e(state)'"
		loc controls = "`e(controls)'"
		loc transvars `free' `state' `controls'
		loc translogNum: word count `transvars'
	
		tempname beta
		mat `beta' = e(b) 
		
		loc freenum: word count `free'
		loc statenum: word count `state'
		loc totnum:  word count `free' `state' 
		
		loc n = 1 /*regenerate the variables used in the routine in order to fit the values*/
		foreach x of local transvars{
			tempvar var_`n' betavar_`n' 
			qui g `betavar_`n'' = `beta'[1,`n'] * `x'
			qui g `var_`n'' = `x'
			loc fit `fit' -`betavar_`n''
			loc ++n
		}
		forv i = 1/`translogNum'{
			forv j = `i'/`translogNum'{ /* `i' */
				tempvar var_`i'`j' betavar_`i'`j'
				cap g `betavar_`i'`j'' = `beta'[1,`n'] * (`var_`i'' * `var_`j'')
				cap g `var_`i'`j'' = (`var_`i'' * `var_`j'')
				loc ++n
			}
		}
		
		/* START - DEFINE BORDERED HESSIAN */
		loc freenum: word count `free'
		loc statenum: word count `state'
		loc totnum:  word count `free' `state' 
		scalar  totnumm = `totnum'
		scalar  totnum_plus_1 = `totnum' + 1
		matrix temp_matrix = J(totnum_plus_1, totnum_plus_1, 0)
		/*insert first derivatives --> first column & first row*/
		foreach a of numlist 1/`totnum' { 
			foreach b of numlist 1/`totnum' { 
								if `a' != `b'	{ 
								cap confirm variable `betavar_`a'`b''
										if !_rc		{
										loc remainder `remainder' + (`betavar_`a'`b''/`var_`a'')
													}
												}
											}
								tempvar fd_`a'
								tempvar first_derivative_`a'
								qui gen `first_derivative_`a'' = (`beta'[1,`a'] + (`betavar_`a'`a''/`var_`a'') `remainder' ) * `e(depvar)' / `var_`a''
								qui su `first_derivative_`a'', meanonly
								loc fd_`a': di %6.3f `r(mean)'
								matrix temp_matrix[1, `= `a' + 1' ] = `fd_`a''
								matrix temp_matrix[`= `a' + 1' , 1] = `fd_`a''
										}
		/*insert second derivatives --> rest of the matrix*/
		foreach a of numlist 1/`totnum' { 
			foreach b of numlist `a'/`totnum' { 
								if `a' != `b'	{ 
								cap confirm variable `betavar_`a'`b''
											if !_rc{
											loc remainder `remainder' + (`betavar_`a'`b'' / `var_`a'')
													}
												}
								tempvar output_elasticity_`b'
								qui gen `output_elasticity_`b'' = (`beta'[1,`b'] + (`betavar_`b'`b'' / `var_`b'') `remainder' )
								tempvar second_derivative_`a'_`b'
								tempvar sd_`a'_`b'
								qui gen `second_derivative_`a'_`b'' = ((`betavar_`a'`b'' / `var_`a'`b'')  + (`output_elasticity_`a'' * `output_elasticity_`b'' ) )   * `e(depvar)' / (`var_`a'' * `var_`b'')							
								qui su `second_derivative_`a'_`b'', meanonly
								loc sd_`a'_`b': di %6.3f `r(mean)'
								matrix temp_matrix[`= `b' + 1' , `= `a' + 1'] = `sd_`a'_`b''
								matrix temp_matrix[`= `a' + 1' , `= `b' + 1'] = `sd_`a'_`b''
										}
									}
				mat list temp_matrix	
				qui gen detbeta = det(temp_matrix) 											
		/*calculate nominator for esa --> first derivatives * quantities */
		matrix mp_vector = J(totnumm ,1, 0) 
		foreach a of numlist 1/`totnum' { 
			foreach b of numlist 1/`totnum' { 
								if `a' != `b'	{ 
								cap confirm variable `betavar_`a'`b''
										if !_rc		{
										loc remainder `remainder' + (`betavar_`a'`b''/`var_`a'')
													}
												}
											}
								tempvar marginal_product_`a' tempvar mp_`a'
								qui gen `marginal_product_`a'' = (`beta'[1,`a'] + (`betavar_`a'`a''/`var_`a'') `remainder' ) * `e(depvar)' / `var_`a''
								qui su `marginal_product_`a'', meanonly
								loc mp_`a': di %6.3f `r(mean)'
								matrix mp_vector[`a',1] = `mp_`a''
									}
		tempname firstde
		mat list  mp_vector
		mata : st_matrix("firstde", colsum(st_matrix("mp_vector")))
		mat list  firstde
		gen fds = .
		replace fds = firstde[1,1] 
		di fds
		/* END - DEFINE BORDERED HESSIAN */
	
		if "`exponential'" != "" {
			qui g `varlist' = exp(`e(depvar)' `fit') `if' // here generate the predicted residuals -- exponential
		}
		else if "`residuals'" != ""{
			qui g `varlist' = `e(depvar)' `fit' `if' // here generate the predicted residuals
		}
		else{ /* in case of 'parameters' option */
			loc freenum: word count `free'
			loc statenum: word count `state'
			loc totnum:  word count `free' `state' 
			forv i = 1/`totnum'{
				forv j = 1/`totnum'{
						if `i' != `j'{ /* generate the cross variables part only  */
						cap confirm variable `betavar_`i'`j''
								if !_rc{
										loc remainder `remainder' + (`betavar_`i'`j''/`var_`i'')
										}
									}
				tempvar betafit_`i'
				qui gen `betafit_`i'' = `beta'[1,`i'] + 2*(`betavar_`i'`i''/`var_`i'') `remainder' // here we use the previously generated variables and weight them by the ith variable
				qui su `betafit_`i'', meanonly
				loc beta_`i': di %6.3f `r(mean)'
									
				/* START - CALCULATE COFACTOR of BORDERED HESSIAN */
				tempvar submatrix
				matrix submatrix = temp_matrix

						forv z = 1/`= `totnum' + 1'{
						matrix submatrix[`z', `i'] = .
						}
								tempname F rF nmF 
								matrix `F' = submatrix
								mata: st_matrix(st_local("rF"), colsum(st_matrix("`F'"))) 
								mata: work = st_matrix("`F'") 
								mata: st_matrix(st_local("nmF"), select(work, colmissing(work) :== 0)) 	
						forv z = 1/`= `totnum' + 1' {
						cap matrix `nmF'[`j', `z'] = .
						}
								tempname cF nnmF
								mata: st_matrix(st_local("cF"), rowsum(st_matrix("`nmF'"))) 
								mata: work = st_matrix("`nmF'") 
								mata: st_matrix(st_local("nnmF"), select(work, rowmissing(work) :== 0)) 
								matrix submatrix = `nnmF'	

				tempvar minor_`i'`j'
				qui gen minor_`i'`j' = det(submatrix)
				
				tempvar cofactor_`i'`j'
				qui gen cofactor_`i'`j' = minor_`i'`j' * ((-1)^(`i'+`j'))  
				/* END - CALCULATE COFACTOR of BORDERED HESSIAN */
				
		tempvar  sigmafit_`i'`j'
		qui gen `sigmafit_`i'`j'' =  (fds * cofactor_`i'`j') / ( `var_`i'' * `var_`j'' * detbeta)
		qui su `sigmafit_`i'`j'', meanonly
		loc sigma_`i'_`j': di %6.3f `r(mean)'
		loc remainder ""
					}
				}
		di _n _n
		di as text "{hline 75}"
		di as text "Translog elasticity estimates" _continue
		di _col(49) "prodest postestimation"
		di as text "{hline 75}"
		di as text "Elasticity Parameter" _continue
		di _col(49) "Value"
		di as text "{hline 75}"
		loc k = 1
		foreach var of varlist `free' `state'{
			di as text "beta_`var'" _continue
			di _col(49) "`beta_`k''"
			loc ++k
		}
		forv i = 1/`translogNum'{
		local x = `i' + 1
				forv j = `x' /`translogNum'{ 
						di as text "sigma_`i'_`j'" _continue
						di _col(49) "`sigma_`i'_`j''"
						}
				}
		di as text "{hline 75}"
		}
	}
end
