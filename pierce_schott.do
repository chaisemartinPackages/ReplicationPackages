use "https://github.com/chaisemartinPackages/ApplicationData/raw/main/pierce_schott_didtextbook.dta", clear

*TWFE regressions
reg delta2001 ntrgap, vce(hc2, dfadjust)
reg delta2002 ntrgap, vce(hc2, dfadjust)
reg delta2004 ntrgap, vce(hc2, dfadjust)
reg delta2005 ntrgap, vce(hc2, dfadjust)

*Test that the NTR-gap treatment is as good as randomly assigned
reg ntrgap lemp1997 lemp1998 lemp1999 lemp2000, vce(hc2, dfadjust)

*Weights analysis
twowayfeweights delta2001 indusid cons ntrgap ntrgap, type(fdTR)

*Heteroscedasticity-robust Yatchew test
yatchew_test delta2001 ntrgap, het_robust
yatchew_test delta2002 ntrgap, het_robust
yatchew_test delta2004 ntrgap, het_robust
yatchew_test delta2005 ntrgap, het_robust

*Parametric linearity test
reg delta2001 ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub
reg delta2002 ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub
reg delta2004 ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub
reg delta2005 ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub

*Pre-trends test
reg delta1999 ntrgap, vce(hc2, dfadjust)
reg delta1998 ntrgap, vce(hc2, dfadjust)
reg delta1997 ntrgap, vce(hc2, dfadjust)

*Pre-trends test, with industry-specific linear trends
reg delta1998lintrend ntrgap, vce(hc2, dfadjust)
reg delta1997lintrend ntrgap, vce(hc2, dfadjust)

*Heteroscedasticity-robust Yatchew test, linear trends
yatchew_test delta2001lintrend ntrgap, het_robust
yatchew_test delta2002lintrend ntrgap, het_robust
yatchew_test delta2004lintrend ntrgap, het_robust
yatchew_test delta2005lintrend ntrgap, het_robust

*Parametric linearity test, linear trends
reg delta2001lintrend ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub
reg delta2002lintrend ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub
reg delta2004lintrend ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub
reg delta2005lintrend ntrgap ntrgapsq ntrgapcub, vce(hc2, dfadjust)
test ntrgapsq ntrgapcub

*Estimators with linear trends
reg delta2001lintrend ntrgap, vce(hc2, dfadjust)
reg delta2002lintrend ntrgap, vce(hc2, dfadjust)
reg delta2004lintrend ntrgap, vce(hc2, dfadjust)
reg delta2005lintrend ntrgap, vce(hc2, dfadjust)