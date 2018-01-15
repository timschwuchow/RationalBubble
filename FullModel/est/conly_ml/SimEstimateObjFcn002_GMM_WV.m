%  Copyright 2012 Timothy John Schwuchow
%  SimEstimateObjFcn002_GMM_WV.m		- Objective function for estimated speculation model ( corresponds to parameterization in brstheory_eq001.pdf)
%  Workign version
function obj	= 	SimEstimateObjFcn002_GMM_WV(par,parset,price,inv,invexp,x,xvar,H,N,b,parest,allpar,moms)


parvec			=	zeros(1,numel(allpar));
ct					=	1;
for i=1:numel(allpar)
	k	=	strmatch(allpar{i},parest);
	if isempty(k)
		parvec(i)	=	parset(i);
	else
		parvec(i)	=	par(ct);
		ct				=	ct+1;
	end
end
%  Unpack
ve		=	parvec(1);
vg		=	parvec(2);
re		=	parvec(3);
rg		=	parvec(4);


%  Generate constants
[vd vi ce1 cg1 cd1 ci1 ce2 cg2 cd2 ci2 cb2 bmvar varv]		=	bayes002([parvec x xvar],b);

[phat dprice ddprice]		=	pricetoerror002(price,inv,invexp,b,x,H,N,bmvar,re,cd2);



T				=	numel(phat);


psi0 			=	ce1 + 1/(1-b*re);
psi1			=	cd1*ce2 - ce1*cd2 - cd2/(1-b*re) - re*ce1;
psi2			=	re*(ce1*cd2 - cd1*ce2);




MomCoefVec 	= 	[psi0 psi1 psi2 zeros(1,moms)];


TheoMoms = [0];
for i=0:moms
	TheoMoms = [TheoMoms sum(ve*MomCoefVec(1:3).*MomCoefVec(1+i:3+i))];
end

ObsMoms = [mean(ddprice)];

for i=0:moms
	cf			=	cov(ddprice(1:T-i-2),ddprice(1+i:T-2));
	ObsMoms		=	[ObsMoms cf(2)];
end

obj			=	mean((TheoMoms-ObsMoms).^2);




