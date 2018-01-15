%  Copyright 2012-2013 Timothy John Schwuchow
%  SimEstimateML_Cons002.m		-	Estimate speculation model using simulated data (maximum likelihood version)
%  Production Version
clear all
if matlabpool('size')
	matlabpool close force local
end
format short g
warning off
path(path,'../../functions');
allpar	=	{'ve','vg','re','rg'};
parest	=	{'ve','vg','re','rg'};
%  Estimation parameters
vegss		=	0.50;
vggss		=	1.00;
regss		=	0.98;
rggss		=	0.500;
pargss	=	[vegss vggss regss rggss];
parub		=	[10.0 10.0 0.999 0.999];
parlb		=	[0.01 0.01 0.100 0.200];
nuse1		=	250;
nuse2		=	1000;
%  Other parameters


NumEst	=	csvread('../../sim/conly/nsprobnum');
NumEst	=	NumEst(2);
matlabpool(8)
for j=1:NumEst
	load(['../../sim/conly/outmat/pricesimconly' int2str(j) '.mat']);
	nuseact1	=	min(nuse1,T);
	nuseact2	=	min(nuse2,T);
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
	smlest		=	zeros(nsim,numel(allpar));
	objsml		=	zeros(nsim,1);
	objvsm		=	zeros(nsim,1);
	objmax		=	zeros(nsim,1);
	optsch		=	optimset('Display','off','TolFun',1e-7,'TolX',1e-10,'MaxIter',30);

	optcon	=	optimset('Display', 'iter','TolFun',1e-7,'TolX',1e-10,'FinDiffType','central','Algorithm','interior-point','UseParallel','always','TolCon',1e-12,'TypicalX',[2 5 0.98 0.7]);

	for i=1:nsim

		tin										=	tic;
		parset									=	[vegrid(j) vggrid(j) regrid(j) rggrid(j)];
		gss										=	parset(parind);
		maxest(i,:)								=	parset;
		smlest(i,:)								=	parset;
		vsmest(i,:)								=	parset;
		[vsmest(i,parind) objvsm(i)]		=	fmincon( @(gss) SimEstimateObjFcn002(gss,parset,pcons(i,1:nuseact1),zeros(1,nuseact1),zeros(1,nuseact1),xgrid(j),xvargrid(j),H,N,b,parest,allpar),gss,[],[],[],[],lb,ub,[],optcon);
		gss										=	vsmest(i,parind);
		if i > 1
			vsmstring 	=	sprintf('True Parameter Value & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nML Estimate & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nDifference & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nStd. Error & (%5.3f) & (%5.3f) & (%5.3f) & (%5.3f) \\\\ \\hline \nTrials & %d &  & &  \\\\ \nObs. per Trial & %d & & & ',parset(1:4),mean(vsmest(1:i,1:4)),parset(1:4) - mean(vsmest(1:i,1:4)),var(vsmest(1:i,1:4)).^0.5,i,nuseact1);
			fprintf('Trial=%d\nObs=%d\n%s\n',i,nuseact1,vsmstring);
%  			vsmfile         =       fopen('tables/VSmallMCEstimates_Cons.tex','w');
%  			fwrite(vsmfile,vsmstring);
%  			fclose(vsmfile);
		end
		[smlest(i,parind) objsml(i)]		=	fmincon( @(gss) SimEstimateObjFcn002(gss,parset,pcons(i,1:nuseact2),zeros(1,nuseact2),zeros(1,nuseact2),xgrid(j),xvargrid(j),H,N,b,parest,allpar),gss,[],[],[],[],lb,ub,[],optcon);
		gss										=	smlest(i,parind);
		if i > 1
			smlstring 	=	sprintf('True Parameter Value & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nML Estimate & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nDifference & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nStd. Error & (%5.3f) & (%5.3f) & (%5.3f) & (%5.3f) \\\\ \\hline \nTrials & %d &  & &  \\\\ \nObs. per Trial & %d & & & ',parset(1:4),mean(smlest(1:i,1:4)),parset(1:4) - mean(smlest(1:i,1:4)),var(smlest(1:i,1:4)).^0.5,i,nuseact2);
			fprintf('Trial=%d\nObs=%d\n%s\n',i,nuseact2,smlstring);
%  			smlfile         =       fopen('tables/SmallMCEstimates_Cons.tex','w');
%  			fwrite(smlfile,smlstring);
%  			fclose(smlfile);
		end
		[maxest(i,parind) objmax(i)]		=	fmincon( @(gss) SimEstimateObjFcn002(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),gss,[],[],[],[],lb,ub,[],optcon);
		if i > 1
			maxstring 	=	sprintf('True Parameter Value & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nML Estimate & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nDifference & %5.3f & %5.3f & %5.3f & %5.3f \\\\ \nStd. Error & (%5.3f) & (%5.3f) & (%5.3f) & (%5.3f) \\\\ \\hline \nTrials & %d &  & &  \\\\ \nObs. per Trial & %d & & &\n',parset(1:4),mean(maxest(1:i,1:4)),parset(1:4) - mean(maxest(1:i,1:4)),var(maxest(1:i,1:4)).^0.5,i,T);
			fprintf('Trial=%d\nObs=%d\n%s\n',i,T,maxstring);
%  			maxfile         =       fopen('tables/LargeMCEstimates_Cons.tex','w');
%  			fwrite(maxfile,maxstring);
%  			fclose(maxfile);
		end
		tout										=	toc(tin);
		fprintf('Iteration %d finished in %5.4f minutes \n',i,tout/60);
%  		if i>1
%  			fprintf('ML Simulation %d :: Trial %d :: Estimate = %5.4f %5.4f %5.4f %5.4f | Obj = %5.4f \n\nTrue                      %5.4f  %5.4f  %5.4f  %5.4f | Obj = %5.4f \nMean                    %5.4f  %5.4f %5.4f %5.4f \nEst-True                   %5.4f %5.4f %5.4f %5.4f\nSt. Err                 (%5.4f) (%5.4f) (%5.4f)  (%5.4f)\nSolved in %5.2f minutes\n****************************\n\n',j,i,maxest(i,1:4), objmax(i),parset,SimEstimateObjFcn002(parset(parind),parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),mean(maxest(1:i,1:4)),mean(maxest(1:i,1:4))-parset(1:4),var(maxest(1:i,1:4)).^0.5,tout/60);
		if i==1
			fprintf('ML Simulation %d :: Trial %d :: Estimate = %5.4f %5.4f %5.4f %5.4f | Obj = %5.4f\nTrue = %5.4f %5.4f %5.4f %5.4f  | Obj = %5.4f\nSolved in %5.4f minutes\n***************\n',j,i,maxest(i,:), objmax(i),parset,SimEstimateObjFcn002(parset(parind),parset,pcons(i,:),zeros(1,T),zeros(1,T),xgrid(j),xvargrid(j),H,N,b,parest,allpar),tout/60);
		end
	end
end
matlabpool close



save('outmat/SimEstimate_Cons	.mat','-V7.3')
quit
