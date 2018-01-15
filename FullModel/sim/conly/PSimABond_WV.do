// Copyright 2013 Timothy John Schwuchow
// PSimABond_WV.m		-	Aralleno-Bond estimator for simulated data
// Working Version


clear all
capture log close _all
set more off
set mem 12G
set matsize 800
timer on 1
local FileName = "PSimABond_WV"
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
	qui insheet using `datdir'PriceSim`x'.csv, comma clear
	qui ren v1 SimNumber
	qui ren v2 t
	qui ren v3 PriceCons
	drop v4
	qui tsset SimNumber t
	forvalues y=1/5	{
		qui xtabond2 D.PriceCons D.L(`y').PriceCons, gmm(D.L(`y').PriceCons, lag(1 3))
		qui matrix b = e(b)
		local est`x'`y'	=	string(b[1,1],"%05.4f")
		local outconly`x' `outconly`x'' & `est`x'`y'
		di "Estimate of lag `y' in simulation `x': " b[1,1]
	}
	!echo "`outconly`x'" >> `logdir'LargeCOnlyResults.tex
}
qui timer off 1
qui timer list 1
di "Program finished running in `r(t1)' seconds"


qui log close `FileName'
