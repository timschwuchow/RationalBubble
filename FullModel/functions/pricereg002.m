%  Copyright 2012 Timothy John Schwuchow
%  pricereg002.m	-	Get regression/correlation coefficients on differenced prices

function [dpreg dpcor]	=	pricereg002(p,lags,pgrid)

regrid	=	pgrid(:,3);
T			=	size(p,2);
nsim		=	size(p,1);
nlags		=	numel(lags);
mlags		=	max(lags);


dpreg		=	zeros(nsim,nlags+1);
dpcor		=	zeros(nsim,nlags);
for i=1:nsim
	pdif					=	(p(i,2:T) - regrid(i)*p(i,1:T-1))';
	TT						=	T-1;
	dpy					=	pdif(1:TT-mlags);
	dpx					=	[ ];
	for t=lags
		dpx				=	[dpx pdif(1+t:TT-mlags+t)];
	end
	dpx					=	[dpx ones(size(dpy))];
	dpreg(i,:)			=	((dpx'*dpx)^(-1)*dpx'*dpy)';
	for t=1:nlags
		cc					=	corrcoef(pdif(1:TT-lags(t)),pdif(1+lags(t):TT));
		dpcor(i,t)		=	cc(2);
	end
end