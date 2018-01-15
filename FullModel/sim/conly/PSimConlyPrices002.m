%  Copyright 2012 Timothy John Schwuchow
%  PSimConlyPrices002.m		-	Solve VF and simulate (consumer-only, corresponds to brstheory_eq001.pdf)
%  Production Version

clear all
format short g
path(path,'../../functions');

load('outmat/specparam.mat');

randn('state',seed);
simstate		=	repng;
x				=	xgrid(repng);
xvar			=	xvargrid(repng);
regrid		=	regrid(repng);
bmvargrid	=	bmvargrid(repng);
varvgrid		=	varvgrid(repng);
cd1grid		=	cd1grid(repng);
ce1grid		=	ce1grid(repng);
vegrid		=	vegrid(repng);
cd2grid		=	cd2grid(repng);
ce2grid		=	ce2grid(repng);
cb2grid		=	cb2grid(repng);






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

save(['outmat/pricesimconly' int2str(repng) '.mat'],'-V7.3');

datout	=	[repng*ones(T,1) (1:T)' pcons(1,:)' peff(1,:)'];
csvwrite(['outcsv/PriceSim' int2str(repng) '.csv'],datout);


%  quit
