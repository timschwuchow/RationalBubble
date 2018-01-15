%  Copyright 2012 Timothy John Schwuchow
%  pricesim002.m	-	Simulate sequences of prices

function [p eta delta inno]	=	pricesim002(pgrid,b,T,H,N,x)

vegrid	=	pgrid(:,1);
vggrid	=	pgrid(:,2);
regrid	=	pgrid(:,3);
rggrid	=	pgrid(:,4);
nsim		=	numel(vegrid);

[vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid]		=	bayes002(pgrid,b);
inno																																		=	innogen002(pgrid,T);

eta					=	zeros(nsim,T);
p						=	zeros(nsim,T);
delta					=	zeros(nsim,T);
eta(:,1)				=	inno(:,1);
for t=1:T
	p(:,t)			=	price002(0,0,eta(:,t),delta(:,t),inno(:,t),x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
	if t < T
		eta(:,t+1)		=	regrid.*eta(:,t) + inno(:,t+1);
		delta(:,t+1)	=	cd2grid.*delta(:,t) + ce2grid.*inno(:,t);
	end
end
