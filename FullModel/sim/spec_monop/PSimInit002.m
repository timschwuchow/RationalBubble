%  Copyright 2012 Timothy John Schwuchow
%  PSimInit002.m		-	Simulate price sequence with speculators (initialization) - t3 includes additional fixed idiosyncratic valuation variance to increase market response to speculator buying
clear all
format short g
path(path,'../../functions');
%  Control parameters
seed		 	=	8;
T				=	1000;
nsim			=	100;
H				=	0.20;
N				=	1;
epts			=	91;
dpts			=	71;
upts			=	55;
ipts			=	60;
esd			=	4.00;
usd			=	3.00;
dsd			=	3.00;
maxshare		=	0.90;
vftolcrit	=	1e-6;
vfmaxiter	=	1000;
howimp		=	25;
howlag		=	2500;
nonaive		=	0;

%  Simulation parameters

x		=	[5.00 15.00 25.00];
xvar	=	[0.000];
ve		=	[1.000 2.000 4.000 8.000];
vg		=	[2.000 5.000 10.00 15.00];
re		=	[0.99];
rg		=	[0.250 0.500 0.750];

%  x		=	[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.500 0.000];
%  xvar	=	[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.500];
%  ve		=	[0.500 0.500 0.500 0.500 0.500 0.500 0.500 1.000 0.250 0.500 0.500];
%  vg		=	[1.000 1.000 1.000 1.000 1.000 2.000 0.500 1.000 1.000 1.000 1.000];
%  re		=	[0.985 0.985 0.985 0.990 0.980 0.985 0.985 0.985 0.985 0.985 0.985];
%  rg		=	[0.700 0.900 0.500 0.700 0.700 0.700 0.700 0.700 0.700 0.700 0.700];
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
[imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid,vdgrid,cb2grid,cd2grid,bmvargrid,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

nsprobout	=	fopen('nsprobnum','w');
fprintf(nsprobout,'%d,%d',ipts,ngrid);
fclose(nsprobout);
save('outmat/specparam.mat','-V7.3');

outfile = fopen('outcsv/parameters.csv','w');
fprintf(outfile,'sim \tve \tvg \tre \trg \tx \tvx \tvd \tvi \tce1 \tcd1 \tce2 \tcd2 \tcb2 \n');
for i=1:ngrid
	fprintf(outfile,['%d' repmat('\t%4.3f',[1 size([pgrid vdgrid vigrid ce1grid cd1grid ce2grid cd2grid cb2grid],2)]) '\n'],[i pgrid(i,:) vdgrid(i) vigrid(i) ce1grid(i) cd1grid(i) ce2grid(i) cd2grid(i) cb2grid(i) bmvargrid(i)]);
end
fclose(outfile);
fprintf('Number of states = %d\n',nstate);

%  quit
