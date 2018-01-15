%  Copyright 2012 Timothy John Schwuchow
%  SimEstimateML_Spec002.m		-	Estimate speculation model using simulated data for speculative model (maximum likelihood version)
%  Production Version
clear all
format short g
%  warning off
path(path,'../../functions');
allpar	=	{'ve','vg','re','rg'};
parest	=	{'ve','vg','re','rg'};
%  Estimation parameters
vegss		=	5.00;
vggss		=	12.0;
regss		=	0.98;
rggss		=	0.6;
pargss	=	[vegss vggss regss rggss];
parub		=	[20.0 20.0 0.999 0.999];
parlb		=	[0.50 0.50 0.700 0.500];
%  Other parameters



load('../../sim/spec/outmat/specparam.mat')
for j=1:ngrid
	load(['../../sim/spec/outmat/pricesim' int2str(j) '.mat']);
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
	mcest		=	zeros(nsim,numel(allpar));
	obj		=	zeros(nsim,1);
	optsch		=	optimset('Display','off','TolFun',1e-6,'TolX',1e-10,'MaxIter',50);
	optcon		=	optimset('Display', 'off','TolFun',1e-6,'TolX',1e-10,'FinDiffType','forward','Algorithm','interior-point','UseParallel','always');
	for i=1:nsim
		parset	=	[vegrid(j) vggrid(j) regrid(j) rggrid(j)];
		mcest(i,:)				=	parset;
		[mcest(i,parind) obj(i)]	=	fminsearch( @(gss) SimEstimateObjFcn002(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),gss,optsch);
		maxest(i,parind)				=	max(min(maxest(i,parind),ub),lb);
		[mcest(i,parind) obj(i)]	=	fmincon( @(gss) SimEstimateObjFcn002(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),gss,[],[],[],[],lb,ub,[],optcon);
		if i>1
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n \nTrue                        %5.4f  %5.4f  %5.4f  %5.4f  \nMean                    %5.4f  %5.4f %5.4f %5.4f \nEst-True                   (%5.4f) (%5.4f) (%5.4f) (%5.4f)\nSt. Dev                 (%5.4f) (%5.4f) (%5.4f)  (%5.4f)\n\n',j,i,mcest(i,1:4), obj(i),parset(1:4),mean(mcest(1:i,1:4)),mean(mcest(1:i,1:4))-parset(1:4),var(mcest(1:i,1:4)).^0.5);
		else
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n \n ',j,i,mcest(i,1:4), obj(i));
		end
	end
end
outstring 	=	sprintf('Mean Estimate & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nTrue Value & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nError & %5.3f & %5.3f & %5.3f & %5.3f & \\\\ \nStd. Error & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nTrials & %d &  & &  \\\\ \nObs. per Trial & %d & & & ',mean(mcest(:,1:4)),parset(1:4),mean(mcest(:,1:4)) - parset(1:4),var(mcest(:,1:4)).^0.5,nsim,T)
outfile		=	fopen('tables/MCEstimates.tex','w');
fwrite(outfile,outstring);
fclose(outfile);

save('outmat/SimEstimate.mat','-V7.3')
quit
