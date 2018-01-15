%  Copyright 2012 Timothy John Schwuchow
%  psimconlymomcheck.m		-	Check validity of moments in generated data


clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');

load('outmat/specparam.mat');
load(['../../sim/conly_t1/outmat/pricesimconly' int2str(1) '.mat']);

re = regrid;
cd1 = cd1grid;
cd2 = cd2grid;
ce1 = ce1grid;
ce2 = ce2grid;
ve  = vegrid;
bmvar = bmvargrid;
bm	=	norminv(1 - H/N,0,bmvar);


j=1;

phat = (1-b*re)*(pcons - (x + bm)/(1-b));

dp = phat(j,2:T) - re*phat(j,1:T-1);
ddp = dp(j,2:T-1) - cd2*dp(j,1:T-2);
%  ddp = reshape(ddp,[1 numel(ddp)]);

momcoef0 = 1 + (1 - b*re)*ce1;
momcoef1 = (1-b*re)*(cd1*ce2-re*ce1) - cd2*(1+(1-b*re)*ce1);
momcoef2 = (1-b*re)*(cd2-re)*cd1*ce2 - cd2*(1-b*re)*(cd1*ce2-re*ce1);
momcoef3 = -(1-b*re)*(cd2-re)*cd1*cd2*ce2;

ErrCoefVec = [momcoef0 momcoef1 momcoef2 momcoef3];

moms = 7;

MomCoefVec = [ErrCoefVec zeros(1,moms)];
TheoMoms = [0];
for i=0:moms
	TheoMoms = [TheoMoms sum(ve*MomCoefVec(1:4).*MomCoefVec(1+i:4+i))];
end

ObsMoms = [mean(ddp)];

for i=0:moms
	cf				=	cov(ddp(1:T-i-2),ddp(1+i:T-2));
	ObsMoms		=	[ObsMoms cf(2)];
end
[TheoMoms;ObsMoms]
obj		=	mean((TheoMoms-ObsMoms).^2)
Y			=	phat(j,2:T)';
X 			= [ones(T-1,1) phat(j,1:T-1)'];
b = (X'*X)^(-1)*X'*Y
e = Y - X*b;
V = (X'*X)^(-1)*X'*X*(X'*X)^-1*var(e)


