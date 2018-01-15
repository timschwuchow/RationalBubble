%  Copyright 2012 Timothy John Schwuchow
%  PSimVF002_WV.m		-	Solve VF and simulate (parallel)
%  Working version
clear all
format short g
path(path,'../../functions');
seed		 	=	8;
T				=	1000;
nsim			=	100;
H				=	0.20;
N				=	1;
epts			=	71;
dpts			=	71;
upts			=	45;
ipts			=	40;
esd			=	4.00;
usd			=	3.00;
dsd			=	3.00;
maxshare		=	0.90;
vftolcrit	=	1e-4;
vfmaxiter	=	1000;
howimp		=	25;
howlag		=	50;
nonaive		=	0;
b				=	0.99;
vdiag			=	0;
NTrials		=	200;
x				=	zeros(NTrials,1);
xvar			=	zeros(NTrials,1);
ve				=	zeros(NTrials,1);
vg				=	zeros(NTrials,1);
re				=	zeros(NTrials,1);
rg				=	zeros(NTrials,1);
pgrid			=	zeros(NTrials,6);
ce1			=	x;
ce2			=	x;
cd1			=	x;
cd2			=	x;
cg1			=	x;
cg2			=	x;
ci1			=	x;
ci2			=	x;
cb2			=	x;
bmvar			=	x;
x(1)		=	5.00;
xvar(1)	=	0.500;
ve(1)		=	1.000;
vg(1)		=	2.000;
re(1)		=	0.99;
rg(1)		=	0.250;
randn('state',seed);

for n=1:NTrials
	fprintf('ve = %5.3f\nvg = %5.3f\nre = %5.3f\nrg = %5.3f\nx = %5.3f\nxvar = %5.3f\n****************\n',ve(n),vg(n),re(n),rg(n),x(n),xvar(n));
	[ve(n) vg(n) re(n) rg(n) x(n) xvar(n) pgrid(n,:) ngrid]		=	grid002(ve(n),vg(n),re(n),rg(n),x(n),xvar(n),b,vdiag);
	[vd(n) vi(n) ce1(n) cg1(n) cd1(n) ci1(n) ce2(n) cg2(n) cd2(n) ci2(n) cb2(n) bmvar(n) varv]		=	bayes002_WV(pgrid(n,:),b);
	[imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid(n,:),vd(n),cb2(n),cd2(n),bmvar(n),ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);
	if n==1
		Vinit		=	repmat(imap',[1  nstate]).*vfprice002(zeros(size(imap)),zeros(size(imap)),emap(smapinv(:,1)),dmap(smapinv(:,3)),umap(smapinv(:,2)),x(n),b,re(n),H,N,bmvar(n),cd1(n),ce1(n));
	else
		Vinit		=	V;
	end

	cnsprob				=	cell(ipts,1);
	t0							=	tic;
	for j=1:ipts
		t1						=	tic;
		[cnsprob{j}]		=	nsprobpar002(pgrid(n,1),pgrid(n,3),ce2(n),cd2(n),cb2(n),bmvar(n),b,H,N,imap,emap,umap,dmap,smapinv,smap,j,nonaive);
		fprintf('Transition matrix %d solved in %5.4f minutes. \n',j,toc(t1)/60);
	end
	fprintf('\n ************************ \nTransition kernel solved in %5.4f seconds \n\n',toc(t0));
	fprintf('Solving policy function \n********************************\n');
	t0				=	tic;
	[pol V]		=	vfsolpar002(re(n),cd1(n),ce1(n),bmvar(n),b,H,N,x(n),vftolcrit,vfmaxiter,cnsprob,imap,emap,dmap,umap,smap,smapinv,howimp,howlag,Vinit);
	paramvec		=	[pgrid(n,:) vd(n) vi(n) ce1(n) cd1(n) ce2(n) cd2(n) cg1(n) cg2(n) ci1(n) ci2(n) bmvar(n)]
	save(['outmat/pol' int2str(n) '.mat'],'pol','V','paramvec','-V7.3');
	clear cnsprob
	if n < NTrials
		ve(n+1)		=	lognrnd(0,0.01)*ve(n);
		vg(n+1)		=	lognrnd(0,0.01)*vg(n);
		re(n+1)		=	max(min(lognrnd(0,0.01)*re(n),0.99),0.50);
		rg(n+1)		=	max(min(lognrnd(0,0.01)*rg(n),0.99),0.00);
		x(n+1)		=	lognrnd(0,0.01)*x(n);
		xvar(n+1)	=	lognrnd(0,0.01)*xvar(n);
	end

end
%  fprintf('\n**************************\n Policy function solved in %5.4f minutes \n',toc(t0)/60);


%  quit

