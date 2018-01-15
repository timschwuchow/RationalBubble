%  Copyright 2012 Timothy John Schwuchow
%  VFSolveDebug002.m		-	Debug VF Solver
%  Test Program
clear all
format short g
path(path,'../../functions');
seed		 	=	8;
T				=	1000;
nsim			=	100;
H				=	0.10;
N				=	1;
epts			=	61;
dpts			=	61;
upts			=	31;
ipts			=	35;
esd			=	4.00;
usd			=	3.00;
dsd			=	3.00;
maxshare		=	0.90;
vftolcrit	=	1e-4;
vfmaxiter	=	1000;
howimp		=	15;
howlag		=	5;
nonaive		=	0;
b				=	0.99;
vdiag			=	0;


x		=	5.00;
xvar	=	20.000;
ve		=	1.000;
vg		=	2.000;
re		=	0.99;
rg		=	0.250;
randn('state',seed);

fprintf('ve = %5.3f\nvg = %5.3f\nre = %5.3f\nrg = %5.3f\nx = %5.3f\nxvar = %5.3f\n****************\n',ve,vg,re,rg,x,xvar);
[ve vg re rg x xvar pgrid ngrid]		=	grid002(ve,vg,re,rg,x,xvar,b,vdiag);
[vd vi ce1 cg1 cd1 ci1 ce2 cg2 cd2 ci2 cb2 bmvar varv]		=	bayes002_WV(pgrid,b);
[imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid,vd,cb2,cd2,bmvar,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

cnsprob				=	cell(ipts,1);
t0							=	tic;
for j=1:ipts
	t1						=	tic;
	[cnsprob{j}]		=	nsprobpar002(ve,re,ce2,cd2,cb2,bmvar,b,H,N,imap,emap,umap,dmap,smapinv,smap,j,nonaive);
	fprintf('Transition matrix %d solved in %5.4f minutes. \n',j,toc(t1)/60);
end
fprintf('\n ************************ \nTransition kernel solved in %5.4f seconds \n\n',toc(t0));
fprintf('Solving policy function \n********************************\n');
V					=	repmat(imap',[1  nstate]).*vfprice002(zeros(size(imap)),zeros(size(imap)),emap(smapinv(:,1)),dmap(smapinv(:,3)),umap(smapinv(:,2)),x,b,re,H,N,bmvar,cd1,ce1);
pol				=	repmat((1:ipts)',1,nstate);
tol				=	1;
iter				=	0;
Vo					=	V;
EV					=	V;
imapstate		=	repmat(imap',1,nstate);
emapstate		=	repmat(emap(smapinv(:,1)),[ipts 1]);
umapstate		=	repmat(umap(smapinv(:,2)),[ipts 1]);
dmapstate		=	repmat(dmap(smapinv(:,3)),[ipts 1]);
pricegrid		=	vfprice002(imap,zeros(size(imap)),emapstate(1,:),dmapstate(1,:),umapstate(1,:),x,b,re,H,N,bmvar,cd1,ce1);
stateextend		=	reshape(repmat((1:nstate),[ipts 1]),[nstate*ipts 1]);
for i=1:25
	fprintf('Iteration %d\n',i);
	Vo					=	V;
	if i>= howlag
		V 				= HowardImp002(pol,stateextend,V,cnsprob,pricegrid,b,imap,imapstate,howimp);
	end
	EV					=	EVEval002(V,cnsprob);
	[V pol]			=	VFMax002(EV,pricegrid,imap,imapstate,b);

	tol				=	max(max(abs(V - Vo)));
	fprintf('Tol = %5.4g\n',tol);
end
flow				=	zeros(ipts,1);
ev					=	flow;
price				=	flow;
s					=	smap(ceil(epts/2),ceil(3*upts/4),ceil(dpts/2));
for i=1:ipts
	price(i)	=	pricegrid(i,s);
	flow(i)	=	(imap(1) - imap(i))*pricegrid(i,s);
	ev(i)		=	b*V(i,:)*cnsprob{i}(:,s);
end
[price imap']
CSVF			=	flow + ev





