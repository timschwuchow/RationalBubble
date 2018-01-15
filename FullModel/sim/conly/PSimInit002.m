%  Copyright 2012-2013 Timothy John Schwuchow
%  PSimInit002.m		-	Simulate price sequence with consumers only (initialization) -  corresponds to parameterization in brstheory_eq001.pdf
%  Production Version
clear all
format short g
path(path,'../../functions');
%  Control parameters
seed		 	=	8;
T				=	10000;
nsim			=	250;
H				=	0.50;
N				=	1;
%  Simulation parameters
x		=	[0.00];
xvar	=	[0.00];
ve		=	[0.50];
vg		=	[1.00];
re		=	[0.98];
rg		=	[0.500];
b		=	0.99;
vdiag	=	0;
%  Optimization parameters
randn('state',seed);

%  Grid out parameter combinations
[vegrid vggrid regrid rggrid xgrid xvargrid pgrid ngrid]		=	grid002(ve,vg,re,rg,x,xvar,b,vdiag);
clear ve vg re rg
%  Generate constants
[vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid]		=	bayes002(pgrid,b);


%  Construct state grid
%  [imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid,vdgrid,cb2grid,cd2grid,bmvargrid,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

nsprobout	=	fopen('nsprobnum','w');
fprintf(nsprobout,'%d,%d',1,ngrid);
fclose(nsprobout);
save('outmat/specparam.mat','-V7.3');
csvwrite('outcsv/pgrid.csv',pgrid);

%  quit
