********************************************************************************
*                                                                              *
*                          MAIN REPLICATION FILE                               *
*       "Two-way Fixed Effects and Differences-in-Differences Estimators       *
*                         with Several Treatments"                             *
*                 de Chaisemartin, C. and D'Haultfoeuille, X.                  *
*                                                                              *
********************************************************************************

// Before running the do files, the users must complete some paths both in the main and the cleaning do files.
// Beware: as the paths are indicated as locals, they must be run simultaneously with the commands they are used by. 

////////// 1. Initialization
clear
set mem 50m
set matsize 120
set more off

** Load Data
use "https://raw.githubusercontent.com/chaisemartinPackages/ApplicationData/main/Hotz_cleaned.dta", clear

********************************************************************************
********************************************************************************

////////// 3. Creation of necessary variables


set matsize 10000

xtset num_st year, delta(5)

** Creating variables necessary for the analysis

	gen d_scrat0dc=d.scrat0dc
	gen F1_scrat0dc=F1.scrat0dc
	gen F1_d_scrat0dc=F.d_scrat0dc
	gen d_edddc=d.edddc
	gen d_ecvalue_per_estab=d.ecvalue_per_estab
	gen F1_edddc=F1.edddc
	gen F2_edddc=F2.edddc
	gen F1_d_edddc=F.d_edddc
	gen L1_edddc=L.edddc
	gen L2_edddc=L2.edddc
	gen F1_d_ecvalue_per_estab=F.d_ecvalue_per_estab



********************************************************************************
********************************************************************************

////////// 2. Analysis

// Dependent variable is the revenue of family home day cares: ecvalue_per_estab


// Treatment variables are:

* edddc : 	    the min years of schooling required to be director of a center-based care,
* scrat0dc :    the minimum staff-to-child ratio,
* no_edddc :    indicator for the fact that there is no minima in terms of years of
			    // schooling required to be director of a center-based care,
* no_scratd0c : indicator for the fact that there is no minima in terms of 
			     // staff-to-child ratio.


** Replication of table 1
	tab edddc
	tab scrat0dc



////// a. Estimation of \hat{\beta_{fe}}, long regression, no controls

** FLAG *** Wrong confidence interval in the paper  
** value of the coefficient on edddc :  -.5658528
**  [95% conf. interval] = [ -.8598469, -.2718588]
	areg ecvalue_per_estab edddc scrat0dc no_scrat0dc no_edddc year1992 year1997, ///
	absorb(num_st) cluster(num_st)


** Showing that \hat{\beta_{fe}} may not be robust to heterogeneous effects

	* When main treatment is the min schooling for directors
	twowayfeweights ecvalue_per_estab nst year edddc, type(feTR) ///
	other_treatments(scrat0dc no_scrat0dc no_edddc) ///
	test_random_weights(pct_black pct_hisp hh_size ln_m_inc college ln_under5 pct_fh_c pct_f_nwork pct_unemploy pct_whome long_comm pct_rural)

	* When main treatment is min staff-to-child ratio
	twowayfeweights ecvalue_per_estab nst year scrat0dc, type(feTR) ///
	other_treatments(no_scrat0dc edddc no_edddc) ///
	test_random_weights(pct_black pct_hisp hh_size ln_m_inc college ln_under5 pct_fh_c pct_f_nwork pct_unemploy pct_whome long_comm pct_rural) 

	* When main treatment is indicator for no min schooling for directors
	//twowayfeweights ecvalue_per_estab nst year no_scrat0dc, type(feTR) ///
	//other_treatments(scrat0dc edddc no_edddc) ///
	//test_random_weights(pct_black pct_hisp hh_size ln_m_inc college ln_under5 pct_fh_c pct_f_nwork pct_unemploy pct_whome long_comm pct_rural) 

	* When main treatment is no min staff-to-child ratio
	//twowayfeweights ecvalue_per_estab nst year no_edddc, type(feTR) ///
	//other_treatments(scrat0dc no_scrat0dc edddc) ///
	//test_random_weights(pct_black pct_hisp hh_size ln_m_inc college ln_under5 pct_fh_c pct_f_nwork pct_unemploy pct_whome long_comm pct_rural) 


	
********************************************************************************

////// b. Estimation of \hat{\beta_{s}}, short regression, no controls
*** FLAG *** Wrong confidence interval in the paper
** value of the coefficient on edddc :   -.0200005
**  [95% conf. interval] = [ -.1167685, .0767675]
	areg ecvalue_per_estab edddc year1992 year1997, absorb(num_st) cluster(num_st)


**** Short regression - computation of the weights
// Here, weights cannot be computed using the twowayfeweights command,
// so they are to be computed by hand. 

	preserve

** Own weights
sum edddc 
scalar mean_D=r(mean)
sum ecvalue_per_estab 
scalar obs=r(N)
gen P_gt=1/obs 
gen nat_weight= P_gt*edddc/mean_D
areg edddc i.year, absorb(nst)
	predict eps_1, residuals 
	gen eps_1_E_D_gt=eps_1*edddc
	sum eps_1_E_D_gt
	scalar denom_W=r(mean)
	gen W=eps_1*mean_D/denom_W
	gen weight=W*nat_weight

egen total_weight_plus=total(weight) if weight>0&weight!=.
egen total_weight_minus=total(weight) if weight<0

	* 56 cells receive a positive weight, sum of the weights is 1.759.
	sum total_weight_plus

	* 71 receive a negative weight, sum of the weights is -0.759.
	sum total_weight_minus

	
** Contamination weights

	* Weights of the effects of not having a requirement on directors' years of schooling
	gen weight_cont1=W*P_gt*no_edddc/mean_D
	gen W_cont1=weight_cont1/P_gt
	gen other_treat1=no_edddc

	egen total_weight_cont1_plus=total(weight_cont1) if weight_cont1>0&weight_cont1!=.
	egen total_weight_cont1_minus=total(weight_cont1) if weight_cont1<0

		* 5 receive a positive weight, they sum up to 0.008,
		sum total_weight_cont1_plus

		* 21 receive a negative weight, they sum to -0.077.
		sum total_weight_cont1_minus


	* Weights of the effects of increasing by one the staff to child ratio
	gen weight_cont3=W*P_gt*scrat0dc/mean_D
	gen W_cont3=weight_cont3/P_gt
	gen other_treat3=scrat0dc

	egen total_weight_cont3_plus=total(weight_cont3) if weight_cont3>0&weight_cont3!=.
	egen total_weight_cont3_minus=total(weight_cont3) if weight_cont3<0

		* 61 effects receive a positive weight, they sum to 0.030
		sum total_weight_cont3_plus
	
		* 87 receive a negative weight, they sum to -0.022.
		sum total_weight_cont3_minus

		
	* Weights of the effects of not having a requirement on staff to child ratio
	gen weight_cont2=W*P_gt*no_scrat0dc/mean_D
	gen W_cont2=weight_cont2/P_gt
	gen other_treat2=no_scrat0dc

	egen total_weight_cont2_plus=total(weight_cont2) if weight_cont2>0&weight_cont2!=.
	egen total_weight_cont2_minus=total(weight_cont2) if weight_cont2<0

		* none receive a positive weight,
		sum total_weight_cont2_plus

		* 5 receive a negative weight, they sum to -0.035.
		sum total_weight_cont2_minus



** Maximum bias of short regression
	gen abs_W_minus_one=abs(W-1)
	sum abs_W_minus_one if edddc>0
	scalar max_bias_own=r(mean)
forvalue i=1/3{
	gen abs_W_cont`i'=abs(W_cont`i')	
	sum abs_W_cont`i' if other_treat`i'>0
	scalar max_bias_other`i'=r(mean)
}
	scalar max_bias_short=max_bias_own+max_bias_other1+max_bias_other2+max_bias_other3
		
restore



**** Long regression - computation of the weights

preserve

** Own weights
sum edddc 
scalar mean_D=r(mean)
sum ecvalue_per_estab 
scalar obs=r(N)
gen P_gt=1/obs 
gen nat_weight= P_gt*edddc/mean_D
areg edddc i.year no_edddc no_scrat0dc scrat0dc, absorb(nst)
	predict eps_1, residuals 
	gen eps_1_E_D_gt=eps_1*edddc
	sum eps_1_E_D_gt
	scalar denom_W=r(mean)
	gen W=eps_1*mean_D/denom_W
	gen weight=W*nat_weight
	
egen total_weight_plus=total(weight) if weight>0&weight!=.
	egen total_weight_minus=total(weight) if weight<0
	sum total_weight_plus
    sum total_weight_minus

** Contamination weights
local j=1
foreach var of varlist no_edddc no_scrat0dc scrat0dc {
	gen weight_cont`j'=W*P_gt*`var'/mean_D
	gen W_cont`j'=weight_cont`j'/P_gt
	gen other_treat`j'=`var'
	local j=`j'+1
	}

forvalue i=1/3{
egen total_weight_cont`i'_plus=total(weight_cont`i') if weight_cont`i'>0&weight_cont`i'!=.
egen total_weight_cont`i'_minus=total(weight_cont`i') if weight_cont`i'<0
sum total_weight_cont`i'_plus
sum total_weight_cont`i'_minus
}

** Maximum bias of long regression
	gen abs_W_minus_one=abs(W-1)
	sum abs_W_minus_one if edddc>0
	scalar max_bias_own=r(mean)
forvalue i=1/3{
	gen abs_W_cont`i'=abs(W_cont`i')	
	sum abs_W_cont`i' if other_treat`i'>0
	scalar max_bias_other`i'=r(mean)
}
	scalar max_bias_long=max_bias_own+max_bias_other1+max_bias_other2+max_bias_other3
	
restore

di max_bias_short, max_bias_long



********************************************************************************

////// c. Alternative estimator

** Studying the sample

// There are 127 (g,t) treated cells (with a non-zero minimum years of schooling)
	tab edddc if edddc!=0

// There are 5 (g,t) cells corresponding to switchers (i.e. in S1)
gen switcher=d_scrat0dc==0&d_edddc!=0&year!=1987&(scrat0dc==0.25|abs(scrat0dc-0.167)<=0.001)&(L1_edddc==0|L1_edddc==12|L1_edddc==14)

	tab switcher
	tab st year if switcher==1


// 122 lost cells, among which some satisfied:

	* First subgroup: 93 (g,t) cells
					  //cells for which the min years of schooling is different 
					  //from zero and never changes
	gen abs_d_edddc=abs(d_edddc)
	bys num_st: egen change_edddc=max(abs_d_edddc)
	tab change_edddc
	tab edddc if change_edddc==0& edddc!=0


	* Second or third subgroup: 19 (g,t) cells
	** Second subgroup: (g,t) cells for which the treatment in t is not 0, 
						//the treatment in t-1 is the same as in t if t is 
						//not the first time period, and the treatment in t+1 
						//is different from the one in t 
						** FLAG *** Is it correct ?
	tab edddc if year==1987&F1_edddc!=edddc&edddc!=0&change_edddc!=0
	tab edddc if year==1992&d_edddc==0&F1_edddc!=edddc&edddc!=0&change_edddc!=0	
	
	** Third subgroup: (g,t) cells for which the treatment is consecutive 
						//only over two time periods, either t and t-1, 
						//or t and t+1 
						** FLAG *** Is it correct ?
	tab edddc if year==1987&F1_edddc==edddc&F2_edddc!=edddc&edddc!=0&change_edddc!=0
	tab edddc if year==1997&L1_edddc==edddc&L2_edddc!=edddc&edddc!=0&change_edddc!=0
	*** If we sum up all the "Total", we have that 4+1+1+13=19

	
	* Fourth or fifth subgroup: 10 (g,t) cells
	** Fourth subggroup: cells (g,t) which become treated from t-1 to t, 
						 //and other treatments change
	tab edddc if year!=1987&d_scrat0dc!=0&d_edddc!=0&edddc!=0&change_edddc!=0
	
	** Fifth subgroup: cells (g,t) which become treated between t-1 and t, 
					   //and other treatments remain stable, but no control (g,t)
					   //cell have the same treatments in t-1 and keeps them at 
					   //the same level
	bys year: tab L1_edddc scrat0dc if edddc!=0&d_edddc!=0&d_edddc!=.&d_scrat0dc==0
	// this line indicates the potential switchers of interest (the 1s in the table)
	bys year: tab L1_edddc scrat0dc if d_edddc==0&d_scrat0dc==0
	// this line indicates the controls
	// to find back (g,t)s of fifth subgroup, one should check that for some switchers, 
	// there is no control with the same values of L1_edddc and scrat0dc.
	// There are 3 switchers in this case. 
	*** If we sum up, we have 7+3=10 (g,t) cells


	
** Computing Forward DIDM 

did_multiplegt ecvalue_per_estab num_st year edddc, if_first_diff(d_scrat0dc==0) trends_nonparam(scrat0dc) always_trends_nonparam breps(0)
ereturn list
// e(N_switchers_effect) indeed equal to 5, 
// e(N_effect_0)!=19+5 because some controls are counted more than once.

scalar didm=e(effect_0)


** Computing Backward DIDM

gen opposite_year=-year
sort num_st year
bys num_st : gen d_scrat0dc_back=scrat0dc[_n+1]-scrat0dc

did_multiplegt ecvalue_per_estab num_st opposite_year edddc, firstdiff_placebo if_first_diff(d_scrat0dc_back==0) trends_nonparam(scrat0dc) always_trends_nonparam breps(0)
ereturn list

scalar didm_f=e(effect_0)


** Computing the verage of forward DIDM and backward DIDM
di (didm+didm_f)/2


e

********************************************************************************
