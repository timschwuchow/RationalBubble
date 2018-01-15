%  Copyright 2012-2013 Timothy John Schwuchow
%  SimEstimateObjFcn002_WV.m		- Objective function for estimated speculation model ( corresponds to parameterization in brstheory_eq001.pdf) (ml version) - Kalman Filter Version
%  Working Version
function obj	= 	SimEstimationObjFcn002_WV(par,parset,price,inv,invexp,x,xvar,H,N,b,parest,allpar,draws)


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
ve			=	parvec(1);
vg			=	parvec(2);
re			=	parvec(3);
rg			=	parvec(4);

[vd vi ce1 cg1 cd1 ci1 ce2 cg2 cd2 ci2 cb2 bmvar varv]		=	bayes002([parvec x xvar],b);
[phat dprice ddprice]													=	pricetoerror002(price,inv,invexp,b,x,H,N,bmvar,re,cd2);

Psi					=	[1/(1-b*re);cd1;ce1];
F						=	[re 0 0 ; 0 cd2 ce2; 0 0 0];
Q						=	[ve 0 ve; 0 0 0; ve 0 ve];
PTTInit				=	[1/(1-re^2)*ve,ce2*re/(1-re*cd2)*ve,ve;ve*ce2*re/(1-re*cd2) vd 0;ve 0 ve];
A						=	chol(PTTInit,'lower');
ndraws				=	size(draws,2);
XiTTInitDraws		=	[zeros(3,1) A*draws];
obj					=	0;
[KalEst KalVar]	=	Kalman002_WV(phat,XiTTInitDraws,PTTInit,Psi,F,Q);
KalVar				=	repmat(KalVar,[1,size(XiTTInitDraws,2)]);
obj					= -mean(sum(log(KalVar.^(-0.5).*exp(-0.5*KalEst.^2.*KalVar.^(-1)))));
%  end
%  for i=1:ndraws
%
%
%  	if isfinite(objtemp)
%  		obj				=	obj + objtemp;
%  	else
%  		obj				=	1e6;
%  	end
%  end
%  obj					=	mean(obj);
%

