%  Copyright 2012-2013 Elliot Anenberg, Patrick Bayer, James Roberts, and Timothy John Schwuchow
%  PSimInitxxx.m		-	Initializes simulation of consumer-only model - defines common control parameters and parameters that vary across simulations
%  Production Version
clear all
format short g
%path(path,'../../functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Control parameters - common to all simulations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seed		 	=	8;			% RNG Seed
nsim			=	1;		% Number of price sequences to be generated for each simulation
T				=	10000;		% Number of periods in each price sequence
H				=	0.50;		% Measure of housing units available in each period
N				=	1.00;		% Measure of consumers who participate in housing market in each period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Parameters - Vary across simulations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x		=	[0.00];				% Fixed individual flow value (population mean)
xvar	=	[0.00];				% Population variance of fixed flow value
ve		=	[0.50 0.90];		% Variance of mu, the innovations in the common flow value component
vg		=	[1.00 1.50];		% Variance of xi, the innovations in the idiosyncratic flow value component
re		=	[0.98];				% rho_eta, the decay parameter of the common flow value component
rg		=	[0.50 0.75];		% rho_gamma, the decay parameter of the common flow value component
b		=	0.99;				% beta, the common discount rate
vdiag	=	0;					% When vdiag=1, the program treats each column of the stacked [x; xvar; ...; rg] vectors as a simulation trial.  Otherwise, each unique element of the cross-product of these vectors is a simulation trial
q       =   0.2;                % proportion of consumers that are noise traders
vz      =   10000;                % variance of noise traders error term

%%%%%%%%%%%%%%%
%%% Set RNG %%%
%%%%%%%%%%%%%%%

rand('state',seed);
randn('state',seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate simulation parameter grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vegrid vggrid regrid rggrid xgrid xvargrid pgrid ngrid]		=	grid002(ve,vg,re,rg,x,xvar,b,vdiag,q,vz);
clear ve vg re rg
%% vegrid, vggrid etc are each vectors from 1:ngrid, with element i of each vector corresponding to the parameter value assumed in the ith simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute coefficients and parameters derived from Bayesian learning process %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid mmgrid]		=	bayes002(pgrid,b,H,vz,q,T);
%  vdgrid = variance of delta (steady state)
%  vigrid = variance of iota (steady state)
%  ce1grid cg1grid cd1grid ci1grid = coefficients on errors in bid functions
%  cd2grid cg2grid cd2grid cb2grid = coefficients in mapping between initial errors and posterior errors after prices are observed
%  bmvargrid - variance of idiosyncratic bid distribution
%  varvgrid = defunct (variance of idiosyncratic bids when consumers are fully informed




%%%%%%%%%%%%
%%% Save %%%
%%%%%%%%%%%%
% nsprobout	=	fopen('nsprobnum','w');
% fprintf(nsprobout,'%d,%d',1,ngrid);
% fclose(nsprobout);
% save('outmat/specparam.mat','-V7.3');
% csvwrite('outcsv/pgrid.csv',pgrid);
PSimConlyPrices002
