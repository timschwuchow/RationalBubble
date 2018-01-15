%  Copyright 2012 Timothy John Schwuchow
%  PSimConlyPrices002.m		-	Solve VF and simulate (consumer-only, corresponds to brstheory_eq001.pdf)
%  Production Version

clear all
format short g
path(path,'../../functions');

load('outmat/specparam.mat');

randn('state',seed);
simstate		=	1;
x				=	xgrid(1);
xvar			=	xvargrid(1);
regrid		=	regrid(1);
bmvargrid	=	bmvargrid(1);
varvgrid		=	varvgrid(1);
cd1grid		=	cd1grid(1);
ce1grid		=	ce1grid(1);
vegrid		=	vegrid(1);
cd2grid		=	cd2grid(1);
ce2grid		=	ce2grid(1);
cb2grid		=	cb2grid(1);






%  Consumer Prices
pcons						=	zeros(nsim,T);
delta						=	zeros(nsim,T);
eta						=	zeros(nsim,T);
u							=	normrnd(0,vegrid^0.5,[nsim T]);
for t=1:T
	if t > 1
		eta(:,t) 	=	regrid*eta(:,t-1) + u(:,t);
		delta(:,t)	=	cd2grid*delta(:,t-1) + ce2grid*u(:,t-1);				% Delta scripted t=t-1 (t is the delta that determines price in period t)
	else
		eta(:,t) =	u(:,t);
	end
	pcons(:,t)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t),delta(:,t),u(:,t),x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
end

peff						=	zeros(nsim,T);
for t=1:T
	peff(:,t)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t),zeros(nsim,1),u(:,t),x,b,regrid,H,N,bmvargrid,0,0);
end

save(['outmat/pricesimconly' int2str(1) '.mat'],'-V7.3');

datout	=	[1*ones(T,1) (1:T)' pcons(1,:)' peff(1,:)'];
csvwrite(['outcsv/PriceSim' int2str(1) '.csv'],datout);


%  quit
