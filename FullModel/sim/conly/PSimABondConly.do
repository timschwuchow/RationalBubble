// Copyright 2013 Timothy John Schwuchow
// PSimABondConly.m		-	Aralleno-Bond estimator for simulated data
// Production Version


clear all
capture log close _all
set more off
set mem 12G
set matsize 800
timer on 1
local FileName = "PSimABondConly"
local datdir = "outcsv/"
local logdir = "logs/"
ashell cat nsprobnum | sed "s/.*,\(.*\)/\1/g"
local NumSim	=	`r(o1)'

qui log using `logdir'`FileName'.txt, text replace name(`FileName')

qui insheet using outcsv/pgrid.csv, clear comma
qui count
forvalues x=1/`r(N)'	{
	local sde`x'	=	string(v1[`x'],"%05.4f")
	local sdg`x'	=	string(v2[`x'],"%05.4f")
	local re`x'		=	string(v3[`x'],"%05.4f")
	local rg`x'		=	string(v4[`x'],"%05.4f")
}

qui set matafavor speed, perm
!if [ -f `logdir'LargeCOnlyResults.tex ]; then rm `logdir'LargeCOnlyResults.tex; fi
forvalues x=1/`NumSim'	{
	local outconly`x'	`sde`x'' & `sdg`x'' & `re`x'' & `rg`x''
	insheet using `datdir'PriceSim`x'.csv, comma clear
	keep in 1/500
	qui ren v1 SimNumber
	qui ren v2 t
	qui ren v3 PriceCons
	drop v4
	qui tsset SimNumber t
	forvalues y=1/5	{
		qui xtabond2 D.PriceCons D.L(`y').PriceCons, gmm(D.L(`y').PriceCons, lag(1 3))
		qui matrix b = e(b)
		local est`x'`y'	=	string(b[1,1],"%05.4f")
		local outconly`x' `outconly`x'' & `est`x'`y''
		di "Estimate of lag `y' in simulation `x': `est`x'`y'' "
	}
	qui sum t, detail
	local outconly`x' `outconly`x'' &  `r(max)'
	local outconly`x' `outconly`x'' \\
	di "`outconly`x''"
	!echo "`outconly`x''" >> `logdir'LargeCOnlyResults.tex
}
qui timer off 1
qui timer list 1
di "Program finished running in `r(t1)' seconds"


qui log close `FileName'
