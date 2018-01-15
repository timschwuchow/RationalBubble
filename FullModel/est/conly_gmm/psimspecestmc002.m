%  Copyright 2012 Timothy John Schwuchow
%  psimspecestmc002.m		-	Estimate speculation model using simulated data (gmm version)
clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');
allpar	=	{'ve','vg','re','rg'};
parest	=	{'ve','vg','re','rg'};
%  Estimation parameters
vegss		=	10.00;
vggss		=	15.0;
regss		=	0.99;
rggss		=	0.50;
pargss	=	[vegss vggss regss rggss];
parub		=	[40.0 40.0 0.999 0.999];
parlb		=	[0.2 0.5 0.01 0.01];
moms		=	10;
%  Other parameters
x			=	1.0;
xvar		=	2.5;
H			=	0.45;
N			=	1;
b			=	0.99;


NumEst	=	csvread('../../sim/conly_t1/nsprobnum');
NumEst	=	NumEst(2);

for j=1:NumEst
	load(['../../sim/conly_t1/outmat/pricesimconly' int2str(j) '.mat']);

	gss		=	[];
	ub			=	[];
	lb			=	[];
	parind	=	[];
	ReEstInd	=	0;
	reIn		=	regrid;
	for i=1:numel(allpar)
		k	=	strmatch(allpar{i},parest);
		if ~ isempty(k)
%  			if allpar{i}=='re'
%  				ReEstInd	=	1;
%  				parest(k) = [];
%  			else
				gss			=	[gss pargss(i)];
				ub				=	[ub parub(i)];
				lb				=	[lb parlb(i)];
				parind		=	[parind i];
%  			end
		end
	end
	%  Options
	mcest		=	zeros(nsim,4);
	obj		=	zeros(nsim,1);
	optcon	=	optimset('Display', 'off','TolFun',1e-8,'Algorithm','interior-point');
	opt		=	optimset('Display', 'off','TolFun',1e-8,'TolX',1e-12);

	for i=1:nsim
%  		if ReEstInd==1
%  			RegEst 	= 	regress(pcons(i,2:T)',[pcons(i,1:T-1)']);
%  			reIn		=	RegEst(1);
%  		end
		parset	=	[vegrid vggrid(j) regrid rggrid(j)];
		mcest(i,:)				=	parset;
		[mcest(i,parind) obj(i)]	=	fmincon( @(gss) psimspecestobj002(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),x,xvar,H,N,b,parest,allpar,moms),gss,[],[],[],[],lb,ub,[],optcon);
		if i>1
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n \n Mean                    %5.4f  %5.4f %5.4f %5.4f \n Est-True                   (%5.4f)  (%5.4f) (%5.4f) (%5.4f)\n St. Dev                  (%5.4f)  (%5.4f) (%5.4f) (%5.4f)\n\n',j,i,mcest(i,:), obj(i),mean(mcest(1:i,:)),mean(mcest(1:i,:))-parset,var(mcest(1:i,:)).^0.5);
		else
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n \n ',j,i,mcest(i,:), obj(i));
		end
	end
%  	save(['out1/psimest' int2str(j) '.mat']);
end
%  [parest2 obj2]	=	fminsearch( @(par) psimspecestobj002(par,pcons(1,:),zeros(1,T),zeros(1,T),x,xvar,H,N,b),pargss,opt);

%  [parest1 obj1]

quit