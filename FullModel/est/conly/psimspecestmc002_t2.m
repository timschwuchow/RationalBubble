%  Copyright 2012 Timothy John Schwuchow
%  psimspecestmc002_t2.m		-	Estimate speculation model using simulated data

clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');
%  Control parameters
seed		 	=	11;
T				=	2000;
nsim			=	200;
x				=	1.0;
xvar			=	2.5;
H				=	0.45;
N				=	1;
epts			=	66;
dpts			=	66;
upts			=	61;
ipts			=	30;
esd			=	2.50;
usd			=	2.50;
dsd			=	2.50;
maxshare		=	0.99;
vftolcrit	=	1e-4;
vfmaxiter	=	1000;
howimp		=	30;
howlag		=	5;
nonaive		=	1;

%  Simulation parameters
ve		=	[3.00];
vg		=	[20.0];
re		=	[0.99];
rg		=	[0.5];
b		=	0.99;
tg		=	100;
vdiag	=	1;
allpar	=	{'ve','vg','re','rg'};
parest	=	{'ve','vg','rg'};

%  Optimization parameters
randn('state',seed);

%  Grid out parameter combinations
[ve vg re rg par ngrid]		=	grid002(ve,vg,re,rg,b,vdiag);

%  Generate constants
[vd vi ce1 cg1 cd1 ci1 ce2 cg2 cd2 ci2 cb2 bmvar varv]		=	bayes002_t1(par,b,xvar);

%  Construct state
[imap emap umap dmap smap smapinv nstate]	=	statemappar002(par,vd,cb2,cd2,bmvar,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

%  Simulate prices

eta		=	zeros(nsim,T,ngrid);
u			=	zeros(nsim,T,ngrid);
delta		=	zeros(nsim,T,ngrid);
pcons		=	zeros(nsim,T,ngrid);


for i=1:ngrid
	fprintf('Simulating data for trial %d \n\n', i)
	u(:,:,i)			=	normrnd(0,ve(i)^0.5,[nsim,T]);
	eta(:,1,i)		=	u(:,1,i);
	pcons(:,1,i)	=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,1,i),delta(:,1,i),u(:,1,i),x,b,re(i),H,N,bmvar(i),cd1(i),ce1(i));
	for t=2:T
		eta(:,t,i)		=	re(i)*eta(:,t-1,i) + u(:,t,i);
		delta(:,t,i)	=	cd2(i)*delta(:,t-1,i) + ce2(i)*u(:,t-1,i);
		pcons(:,t,i)	=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t,i),delta(:,t,i),u(:,t,i),x,b,re(i),H,N,bmvar(i),cd1(i),ce1(i));
	end
end


%  Estimation Parameters
vegss		=	5;
vggss		=	12;
regss		=	re;
rggss		=	rg;
pargss	=	[vegss vggss regss rggss];
parub		=	[40.0 40.0 0.999 0.999];
parlb		=	[0.2 0.5 0.01 0.01];

gss		=	[];
ub			=	[];
lb			=	[];
parset	=	[];
parind	=	[];
for i=1:numel(allpar)
	k	=	strmatch(allpar{i},parest);
	if ~ isempty(k)
		gss			=	[gss pargss(i)];
		ub				=	[ub parub(i)];
		lb				=	[lb parlb(i)];
		parind		=	[parind i];
	end
end



for i=1:ngrid
	fprintf('Estimating trial %d \n **************************\n',i);
	%  Options
	mcest		=	zeros(nsim,4);
	obj		=	zeros(nsim,1);
	optcon	=	optimset('Display', 'off','Algorithm','interior-point');
	opt		=	optimset('Display', 'off','TolFun',1e-8,'TolX',1e-12);
	parset	=	[ve(i) vg(i) 0.999 rg(i)];
	for j=1:nsim
		mcest(j,:)						=	parset;
		[mcest(j,parind) obj(j)]	=	fmincon( @(gss) psimspecestobj002(gss,parset,pcons(j,:,i),zeros(1,T),zeros(1,T),x,xvar,H,N,b,parest,allpar),gss,[],[],[],[],lb,ub,[],optcon);
		if j>1
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4g \n \n Mean                    %5.4f  %5.4f %5.4f %5.4f \n St. Dev                  (%5.4f)  (%5.4f) (%5.4f) (%5.4f)\n\n',i,j,mcest(j,:), obj(j),mean(mcest(1:j,:)),var(mcest(1:j,:)).^0.5);
		else
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4g \n \n ',i,j,mcest(j,:), obj(j));
		end
	end
end
%  [parest2 obj2]	=	fminsearch( @(par) psimspecestobj002(par,pcons(1,:),zeros(1,T),zeros(1,T),x,xvar,H,N,b),pargss,opt);

%  [parest1 obj1]

