%  Copyright 2012-2013 Timothy John Schwuchow
%  SimEstimateML_Cons002.m		-	Estimate speculation model using simulated data (maximum likelihood version)
%  Working Version
clear all
format short g
warning off
path(path,'../../functions');
allpar	=	{'ve','vg','re','rg'};
parest	=	{'ve','vg','re','rg'};
%  Estimation parameters
vegss		=	6.00;
vggss		=	10.0;
regss		=	0.985;
rggss		=	0.600;
pargss	=	[vegss vggss regss rggss];
parub		=	[25.0 25.0 0.999 0.999];
parlb		=	[0.01 0.01 0.400 0.500];
ndraws	=	250;
seed		=	15;
randn('state',seed);
%  Other parameters


NumEst	=	csvread('../../sim/conly/nsprobnum');
NumEst	=	NumEst(2);
draws		=	normrnd(0,1,[3,ndraws]);
for j=1:NumEst
	load(['../../sim/conly/outmat/pricesimconly' int2str(j) '.mat']);
	gss		=	[];
	ub			=	[];
	lb			=	[];
	parind	=	[];
	reIn		=	regrid;
	for i=1:numel(allpar)
		k	=	strmatch(allpar{i},parest);
		if ~ isempty(k)
			gss			=	[gss pargss(i)];
			ub				=	[ub parub(i)];
			lb				=	[lb parlb(i)];
			parind		=	[parind i];
		end
	end
	%  Options
	maxest		=	zeros(nsim,numel(allpar));
	objmax		=	zeros(nsim,1);
	optsch		=	optimset('Display','iter','TolFun',1e-7,'TolX',1e-10,'MaxIter',20);
%  	optcon		=	optimset('Display', 'off','TolFun',1e-8,'TolX',1e-10,'FinDiffType','central','Algorithm','interior-point','UseParallel','always');

%  	optcon		=	optimset('Display', 'iter','TolFun',1e-10,'TolX',1e-12,'FinDiffType','central','Algorithm','interior-point','UseParallel','always','TolCon',1e-12);
	optcon	=	optimset('Display', 'iter','TolFun',1e-7,'TolX',1e-10,'FinDiffType','central','Algorithm','interior-point','UseParallel','always','TolCon',1e-12);
%  	optsqp	=	optimset('Display', 'off','TolFun',1e-10,'TolX',1e-12,'FinDiffType','central','Algorithm','sqp','UseParallel','always','TolCon',1e-12,'Hessian','bfgs','LargeScale','off');
%  	optga			=	gaoptimset('Display', 'off','UseParallel','always','PopulationSize',100,'Generations',10);

	for i=1:nsim
		tin										=	tic;
		parset									=	[vegrid(j) vggrid(j) regrid(j) rggrid(j)];
		maxest(i,:)								=	parset;
%  		[maxest(i,parind) objmax(i)]		=	ga( @(gss) SimEstimateObjFcn002_WV2(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),numel(parind),[],[],[],[],lb,ub,[],optga);
%  		[maxest(i,parind) objmax(i)]		=	fminsearch( @(gss) SimEstimateObjFcn002(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),gss,optsch);
%  		maxest(i,parind)						=	max(min(maxest(i,parind),ub),lb);
%  		gss										=	maxest(i,parind);

		[maxest(i,parind) objmax(i)]		=	fmincon( @(gss) SimEstimateObjFcn002_WV(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar,draws),gss,[],[],[],[],lb,ub,[],optcon);
		tout										=	toc(tin);
		if i>1
			fprintf('ML Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n\nTrue                      %5.4f  %5.4f  %5.4f  %5.4f \nMean                    %5.4f  %5.4f %5.4f %5.4f \nEst-True                   %5.4f %5.4f %5.4f %5.4f\nSt. Err                 (%5.4f) (%5.4f) (%5.4f)  (%5.4f)\nSolved in %5.2f minutes\n****************************\n\n',j,i,maxest(i,1:4), objmax(i),parset,mean(maxest(1:i,1:4)),mean(maxest(1:i,1:4))-parset(1:4),var(maxest(1:i,1:4)).^0.5,tout/60);
		else
			fprintf('ML Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f\nSolved in %5.4f minutes\n***************\n',j,i,maxest(i,:), objmax(i),tout/60);
		end
	end
end
outstring 	=	sprintf('True Parameter Value & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nML Estimate & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nStd. Error & (%5.3f) & (%5.3f) & (%5.3f) & (%5.3f) \\\\ \\hline \nTrials & %d &  & &  \\\\ \nObs. per Trial & %d & & & ',parset(1:4),mean(maxest(:,1:4)),var(maxest(:,1:4)).^0.5,nsim,T)
outfile		=	fopen('tables/MCEstimates_Cons.tex','w');
fwrite(outfile,outstring);
fclose(outfile);

save('outmat/SimEstimate_Cons	.mat','-V7.3')
quit
