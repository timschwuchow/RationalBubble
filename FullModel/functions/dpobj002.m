%  Copyright 2012 Timothy John Schwuchow
%  dpobj002.m	-	Objective function for moments of dphat

function obj	=	dpobj002(pg,pfix,p,x,b,H,N,lags,bmvar,varsch,varfix,ew)


%  Unpack
nv					=	max(max(varsch),max(varfix));
pvec				=	zeros(1,nv);
pvec(1,varfix)	=	pfix;
pvec(1,varsch)	=	pg;
ve					=	pvec(1);
vg					=	pvec(2);
re					=	pvec(3);
rg					=	pvec(4);

[vd vi ce1 cg1 cd1 ci1 ce2 cg2 cd2 ci2 cb2 bmvar]		=	bayes002(pvec,b);
nlags				=	numel(lags);
mlags				=	max(lags);
T					=	size(p,2);
dphat				=	zeros(1,T-1);

%  Construct dphat
bm				=	norminv(1-H/N,0,bmvar);
p1				=	p - (x + bm)/(1-b);
dphat			=	(1-b*re)*(p1(2:T) - re*p1(1:T-1));

momset		=	zeros(nlags+1,T-1);
momset(1,:)	=	dphat;
for t=1:nlags
	momset(t+1,1:T-1-lags(t))		=	dphat(1:T-1-lags(t)).*dphat(1+lags(t):T-1) - dphatmom002(lags(t),pvec,vd,ce1,ce2,cd1,cd2,b);
end
momsubset	=	momset(:,1:T-1-mlags);
momvect		=	mean(momsubset,2)';
if ew
	obj				=	T*momvect*(momsubset*momsubset')^(-1)*momvect';
else
	obj				=	T*momvect*eye(numel(momvect))*momvect';
end