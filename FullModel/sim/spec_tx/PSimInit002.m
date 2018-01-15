%  Copyright 2012-2013 Timothy John Schwuchow
%  PSimInit002_WV.m		-	Simulate price sequence with speculators (initialization)
%  Working version
clear all
format short g
path(path,'../../functions');
%  Control parameters
seed		 	=	8;
T				=	1000;
nsim			=	100;
H				=	0.40;
N				=	1;
epts			=	85;
dpts			=	85;
upts			=	55;
ipts			=	55;
esd			=	4.00;
usd			=	3.00;
dsd			=	3.00;
maxshare		=	0.90;
vftolcrit	=	1e-6;
vfmaxiter	=	5000;
howimp		=	25;
howlag		=	25;
nonaive		=	0;

%  Simulation parameters

x		=	[5.000 15.00];
xvar	=	[0.000 100.0];
ve		=	[1.000];
vg		=	[1.200 0.800 1.000 1.100 1.400 1.600 1.800 2.000 3.000];
re		=	[0.990 0.980];
rg		=	[0.500 0.900];
tx		=	0.03;

b		=	0.99;
vdiag	=	0;
%  Optimization parameters
randn('state',seed);

%  Grid out parameter combinations
[vegrid vggrid regrid rggrid xgrid xvargrid pgrid ngrid]		=	grid002(ve,vg,re,rg,x,xvar,b,vdiag);
clear ve vg re rg x xvar
%  Generate constants
[vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid]		=	bayes002_WV(pgrid,b);



%  Construct state grid
[imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid,vdgrid,cb2grid,cd2grid,bmvargrid,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

nsprobout	=	fopen('nsprobnum','w');
fprintf(nsprobout,'%d,%d',ipts,ngrid);
fclose(nsprobout);
save('outmat/specparam.mat','-V7.3');

outfile = fopen('outcsv/parameters.csv','w');
fprintf(outfile,'sim \tve \tvg \tre \trg \tx \tvx \tvd \tvi \tce1 \tcd1 \tce2 \tcd2 \tcb2 \tbmvar \n');
for i=1:ngrid
	fprintf(outfile,['%d' repmat('\t%4.3f',[1 size([pgrid vdgrid vigrid ce1grid cd1grid ce2grid cd2grid cb2grid bmvargrid],2)]) '\n'],[i pgrid(i,:) vdgrid(i) vigrid(i) ce1grid(i) cd1grid(i) ce2grid(i) cd2grid(i) cb2grid(i) bmvargrid(i)]);
end
fclose(outfile);
fprintf('Number of states = %d | Number of simulations %d\n',nstate,ngrid);


quit
