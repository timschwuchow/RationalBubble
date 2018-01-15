%  Copyright 2012 Timothy John Schwuchow
%  PSimPrices002.m		-	Simulate Prices
%  Production Version

clear all
format short g
path(path,'../../functions');

load('outmat/specparam.mat');
load(['outmat/pol' int2str(repng) '.mat']);

randn('state',seed);
emap			=	emap(repng,:);
dmap			=	dmap(repng,:);
umap			=	umap(repng,:);
regrid		=	regrid(repng);
bmvargrid	=	bmvargrid(repng);
varvgrid		=	varvgrid(repng);
cd1grid		=	cd1grid(repng);
ce1grid		=	ce1grid(repng);
vegrid		=	vegrid(repng);
cd2grid		=	cd2grid(repng);
ce2grid		=	ce2grid(repng);
cb2grid		=	cb2grid(repng);
x				=	xgrid(repng);

%  Monopolist Speculative Prices

PriceSpecMon			=	zeros(nsim,T);
StateSpecMon			=	zeros(nsim,T);
[~,einit]				=	min(abs(emap - 0));
[~,uinit]				=	min(abs(umap - 0));
[~,dinit]				=	min(abs(dmap - 0));
StateSpecMon(:,1)		=	smap(einit,uinit,dinit);
SpecInvMonState		=	zeros(nsim,T);
SpecInvMon				=	zeros(nsim,T);
iinit						=	1;
cf							=	zeros(nsim,T);
for t=1:T
	if t==1
		SpecInvMonState(:,t)		=	reshape(pol(iinit,StateSpecMon(:,t)),[nsim 1]);
	else
		SpecInvMonState(:,t)		=	reshape(pol(sub2ind(size(pol),SpecInvMonState(:,t-1),StateSpecMon(:,t))),[nsim 1]);
	end
	SpecInvMon(:,t)		=	imap(SpecInvMonState(:,t))';
	PriceSpecMon(:,t)		=	pricespec002(SpecInvMon(:,t),zeros(nsim,1),emap(smapinv(StateSpecMon(:,t),1))',dmap(smapinv(StateSpecMon(:,t),3))',umap(smapinv(StateSpecMon(:,t),2))',x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);

	if t < T
		naiveerr				=	(norminv(1 - (H)/ N,0,bmvargrid) - norminv(1 - (H - SpecInvMon(:,t))/ N,0,bmvargrid));
		[~,utemp]			=	min(abs(repmat(umap,nsim,1) - repmat(normrnd(0,vegrid^0.5,[nsim 1]),[1 upts])),[],2);
		[~,etemp]			=	min(abs(repmat(emap,nsim,1) - (repmat(regrid*emap(smapinv(StateSpecMon(:,t),1))' + umap(utemp)',[1 epts]))),[],2);
		[~,dtemp]			=	min(abs(repmat(dmap,nsim,1) - (repmat(cb2grid*naiveerr + cd2grid*dmap(smapinv(StateSpecMon(:,t),3))' + ce2grid*umap(smapinv(StateSpecMon(:,t),2))',[1 dpts]))),[],2);
		StateSpecMon(:,t+1)	=	smap(sub2ind(size(smap),etemp,utemp,dtemp));
	end
end


% Competitive Speculator Prices
PriceSpecComp			=	zeros(nsim,T);
StateSpecComp			=	zeros(nsim,T);
StateSpecComp(:,1)	=	smap(einit,uinit,dinit);
SpecInvCompState		=	zeros(nsim,T);
SpecInvComp				=	zeros(nsim,T);
thresh					=	1-1e-10;
for t=1:T
	if t > 1
		utemp						=	smapinv(StateSpecMon(:,t),2);
		etemp						=	smapinv(StateSpecMon(:,t),1);
		[~,dtemp]				=	min(abs(repmat(dmap,nsim,1) - (repmat(cd2grid*dmap(smapinv(StateSpecComp(:,t-1),3))' + ce2grid*umap(smapinv(StateSpecComp(:,t-1),2))' + cb2grid*(norminv(1-H/N,0,bmvargrid) - norminv(1-(H-SpecInvComp(:,t-1))/N,0,bmvargrid)),[1 dpts]))),[],2);
		StateSpecComp(:,t)	=	smap(sub2ind(size(smap),etemp,utemp,dtemp));
		SpecInvComp(:,t)		=	max(N*(normcdf( -(x+b*norminv(1-H/N,0,bmvargrid))/(1-b) - emap(smapinv(StateSpecComp(:,t),1))'/(1-b*regrid) - cd1grid*dmap(smapinv(StateSpecComp(:,t),3))' - ce1grid*umap(smapinv(StateSpecComp(:,t),2))',0,bmvargrid)-1)+H*thresh,0);
	else
		SpecInvComp(:,t)			=	max(N*(normcdf( -(x+b*norminv(1-H/N,0,bmvargrid))/(1-b) - emap(smapinv(StateSpecComp(:,t),1))'/(1-b*regrid) - cd1grid*dmap(smapinv(StateSpecComp(:,t),3))' - ce1grid*umap(smapinv(StateSpecComp(:,t),2))',0,bmvargrid)-1)+H*thresh,0);
	end
	PriceSpecComp(:,t)		=	pricespec002(SpecInvComp(:,t),zeros(nsim,1),emap(smapinv(StateSpecComp(:,t),1))',dmap(smapinv(StateSpecComp(:,t),3))',umap(smapinv(StateSpecComp(:,t),2))',x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
end

%  Consumer Prices

StateCons					=	zeros(nsim,T);
StateCons(:,1)				=	StateSpecMon(:,1);
PriceCOnly					=	zeros(nsim,T);
for t=1:T
	PriceCOnly(:,t)		=	pricespec002(zeros(nsim,1),zeros(nsim,1),emap(smapinv(StateCons(:,t),1))',dmap(smapinv(StateCons(:,t),3))',umap(smapinv(StateCons(:,t),2))',x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
	if t < T
		utemp					=	smapinv(StateSpecMon(:,t+1),2);
		etemp					=	smapinv(StateSpecMon(:,t+1),1);
		[~,dtemp]			=	min(abs(repmat(dmap,nsim,1) - (repmat(cd2grid*dmap(smapinv(StateCons(:,t),3))' + ce2grid*umap(smapinv(StateCons(:,t),2))',[1 dpts]))),[],2);
		StateCons(:,t+1)	=	smap(sub2ind(size(smap),etemp,utemp,dtemp));
	end
end

PriceEff						=	zeros(nsim,T);
for t=1:T
	PriceEff(:,t)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),emap(smapinv(StateCons(:,t),1))',zeros(nsim,1),umap(smapinv(StateCons(:,t),2))',x,b,regrid,H,N,bmvargrid,0,0);
end
save(['outmat/pricesim' int2str(repng) '.mat'],'PriceSpecMon','PriceSpecComp','PriceCOnly','PriceEff','SpecInvMon','StateSpecMon','StateCons','StateSpecComp', '-V7.3');

OutPrices				=	[repng*ones(T,1) (1:T)' PriceSpecMon(1,:)' PriceSpecComp(1,:)' PriceCOnly(1,:)' PriceEff(1,:)' SpecInvMon(1,:)' SpecInvComp(1,:)' emap(smapinv(StateCons(1,:),1))' umap(smapinv(StateCons(1,:),2))' dmap(smapinv(StateSpecMon(1,:),3))' dmap(smapinv(StateSpecComp(1,:),3))' dmap(smapinv(StateCons(1,:),3))' repmat(pgrid(repng,:),[T 1])];
csvwrite(['outcsv/Pricerepng.csv'],OutPrices);


