%  Copyright 2013 Timothy John Schwuchow
%  SimPricesSpecComp.m		-	Simulate Prices for competitive speculative market
%  Production Version

clear all
format short g
path(path,'../../functions');

load('outmat/specparam.mat');

randn('state',seed);
for i=1:size(pgrid,1)
	repng			=	i;
	simstate		=	repng;
	x				=	xgrid(repng);

	xvar			=	xvargrid(repng);
	re				=	regrid(repng);
	bmvar			=	bmvargrid(repng);
	cd1			=	cd1grid(repng);
	ce1			=	ce1grid(repng);
	ve				=	vegrid(repng);
	cd2			=	cd2grid(repng);
	ce2			=	ce2grid(repng);
	cb2			=	cb2grid(repng);

	%  bmvargrid	=	150;


	PriceSpecComp					=	zeros(nsim,T);
	PriceCons						=	zeros(nsim,T);
	DeltaSpec						=	zeros(nsim,T);
	DeltaCons						=	zeros(nsim,T);
	Eta								=	zeros(nsim,T);
	SpecInv							=	zeros(nsim,T);
	SpecInvRaw						=	zeros(nsim,T);
	Score								=	zeros(nsim,T);
	Mu									=	normrnd(0,ve^0.5,[nsim T]);

	for t=1:T
		if t > 1
			Eta(:,t) 				=	re*Eta(:,t-1) + Mu(:,t);
			DeltaSpec(:,t)			=	cd2*DeltaSpec(:,t-1) + ce2*Mu(:,t-1);				% Delta scripted t=t-1 (t is the delta that determines price in period t)
			DeltaCons(:,t)			=	cd2*DeltaSpec(:,t-1) + ce2*Mu(:,t-1);				% Delta scripted t=t-1 (t is the delta that determines price in period t)
			SpecInvRaw(:,t)		=	0.99*N*(normcdf( -x/(1-b) -b*norminv(1-H/N,0,bmvar)/(1-b) - Eta(:,t)/(1-b*re) + (b*cd1*ce2-ce1)*Mu(:,t) + cd1*(b*cd2 - 1)/(1-b*cd2)*(DeltaSpec(:,t)+ce2/cd2*Mu(:,t)),0,bmvar)-1)+H*0.99;
			SpecInv(:,t)			=	max(N*(normcdf( -x/(1-b) -b*norminv(1-H/N,0,bmvar)/(1-b) - Eta(:,t)/(1-b*re) + (b*cd1*ce2-ce1)*Mu(:,t) + cd1*(b*cd2 - 1)/(1-b*cd2)*(DeltaSpec(:,t)+ce2/cd2*Mu(:,t)),0,bmvar)-1)+H*0.99,0);
			Score(:,t)				=	-x/(1-b) -b*norminv(1-H/N,0,bmvar)/(1-b) - Eta(:,t)/(1-b*re) + (b*cd1*ce2-ce1)*Mu(:,t) + cd1*(b*cd2 - 1)/(1-b*cd2)*(DeltaSpec(:,t)+ce2/cd2*Mu(:,t));
		else
			Eta(:,t) 				=	Mu(:,t);
			SpecInvRaw(:,t)		=	0.99*N*(normcdf( -x/(1-b) -b*norminv(1-H/N,0,bmvar)/(1-b) - Eta(:,t)/(1-b*re) + (b*cd1*ce2-ce1)*Mu(:,t) + cd1*(b*cd2 - 1)/(1-b*cd2)*(ce2/cd2*Mu(:,t)),0,bmvar)-1)+H*0.99;
			SpecInv(:,t)			=	max(N*(normcdf( -x/(1-b) -b*norminv(1-H/N,0,bmvar)/(1-b) - Eta(:,t)/(1-b*re) + (b*cd1*ce2-ce1)*Mu(:,t) + cd1*(b*cd2 - 1)/(1-b*cd2)*(ce2/cd2*Mu(:,t)),0,bmvar)-1)+H*0.99,0);
			Score(:,t)				=	-x/(1-b) -b*norminv(1-H/N,0,bmvar)/(1-b) - Eta(:,t)/(1-b*re) + (b*cd1*ce2-ce1)*Mu(:,t) + cd1*(b*cd2 - 1)/(1-b*cd2)*(ce2/cd2*Mu(:,t));
		end
		PriceSpecComp(:,t)		=	pricespec002(SpecInv(:,t),zeros(nsim,1),Eta(:,t),DeltaSpec(:,t),Mu(:,t),x,b,re,H,N,bmvar,cd1,ce1);
	end
	PriceDif							=	b*PriceSpecComp(:,2:T) - PriceSpecComp(:,1:T-1);
	PriceInd							=	(SpecInv(:,2:T)>0).*(SpecInv(:,2:T)<0.99*H).*(SpecInv(:,1:T-1)>0).*(SpecInv(:,1:T-1)<0.99*H);

	fprintf('Parameters : %5.3f  %5.3f  %5.3f  %5.3f   %5.3f  %5.3f  \nMean Price Difference %5.3f  %5.3f \nMinInv %d  MaxInv %d\n************************\n',pgrid(i,:),mean(mean(PriceDif(PriceInd==1))),mean(mean(PriceDif)),sum(sum(SpecInv==min(min(SpecInv)))),sum(sum(SpecInv==max(max(SpecInv)))))
end
%  Consumer Prices
%  pcons						=	zeros(nsim,T);
%
%  for t=1:T
%  	if t > 1
%  		eta(:,t) 	=	regrid*eta(:,t-1) + u(:,t);
%  		delta(:,t)	=	cd2grid*delta(:,t-1) + ce2grid*u(:,t-1);				% Delta scripted t=t-1 (t is the delta that determines price in period t)
%  	else
%  		eta(:,t) =	u(:,t);
%  	end
%  	pcons(:,t)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t),delta(:,t),u(:,t),x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
%  end
%
%  peff						=	zeros(nsim,T);
%  for t=1:T
%  	peff(:,t)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t),zeros(nsim,1),u(:,t),x,b,regrid,H,N,bmvargrid,0,0);
%  end
%
%  save(['outmat/pricesimconly' int2str(repng) '.mat'],'-V7.3');
%
%  datout	=	[repng*ones(T,1) (1:T)' pcons(1,:)' peff(1,:)'];
%  csvwrite(['outcsv/PriceSim' int2str(repng) '.csv'],datout);
%

%  quit
