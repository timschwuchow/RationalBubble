%  Copyright 2012 Timothy John Schwuchow
%  PSimInit002.m		-	Simulate price sequence with speculators (initialization) - t3 includes additional fixed idiosyncratic valuation variance to increase market response to speculator buying
clear all
format short g
path(path,'../../functions');
%  Control parameters
seed		 	=	11;
T				=	1000;
nsim			=	500;
H				=	0.50;
N				=	1;
epts			=	85;
dpts			=	85;
upts			=	45;
ipts			=	40;
esd			=	3.5;
usd			=	4;
dsd			=	1.5;
maxshare		=	0.99;
vftolcrit		=	1e-4;
vfmaxiter		=	1000;
howimp			=	45;
howlag			=	2;
nonaive			=	0;

%  Simulation parameters
x		=	[0.00];
xvar	=	[2.00];
ve		=	[0.50 1.00 0.25];
vg		=	[1.00 2.00 0.50];
re		=	[0.99 0.98];
rg		=	[0.70 0.50];
b		=	0.99;
vdiag		=	0;
%  Optimization parameters
randn('state',seed);

%  Grid out parameter combinations
[vegrid vggrid regrid rggrid xgrid xvargrid pgrid ngrid]		=	grid002(ve,vg,re,rg,x,xvar,b,vdiag);
clear ve vg re rg
%  Generate constants
[vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid]		=	bayes002(pgrid,b);



%  Construct state grid
[imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid,vdgrid,cb2grid,cd2grid,bmvargrid,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

nsprobout	=	fopen('nsprobnum','w');
fprintf(nsprobout,'%d,%d',ipts,ngrid);
fclose(nsprobout);
save('outmat/specparam.mat','-V7.3');

csvwrite('outcsv/pgrid.csv',pgrid);

%  quit
