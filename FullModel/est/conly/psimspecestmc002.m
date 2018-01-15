%  Copyright 2012 Timothy John Schwuchow
%  psimspecestmc002.m		-	Estimate speculation model using simulated data
clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');
allpar	=	{'ve','vg','re','rg'};
parest	=	{'ve'};
ngrid		=	5;
%  Estimation parameters
vegss		=	8.00;
vggss		=	15.0;
regss		=	0.98;
rggss		=	0.50;
pargss	=	[vegss vggss regss rggss];
parub		=	[40.0 40.0 0.999 0.999];
parlb		=	[0.2 0.5 0.01 0.01];
moms		=	2;
%  Other parameters
x			=	1.0;
xvar		=	2.5;
H			=	0.45;
N			=	1;
b			=	0.99;




for j=1:1
	load(['../../sim/conly/outmat/pricesimconly' int2str(j) '.mat']);

	gss		=	[];
	ub			=	[];
	lb			=	[];
	parind	=	[];
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
	mcest		=	zeros(nsim,4);
	obj		=	zeros(nsim,1);
	optcon	=	optimset('Display', 'off','TolFun',1e-8,'Algorithm','interior-point');
	opt		=	optimset('Display', 'off','TolFun',1e-8,'TolX',1e-12);
	parset	=	[vegrid vggrid(j) regrid rggrid(j)];
	for i=1:nsim
		mcest(i,:)				=	parset;
		[mcest(i,parind) obj(i)]	=	fmincon( @(gss) psimspecestobj002(gss,parset,pcons(i,:),zeros(1,T),zeros(1,T),x,xvar,H,N,b,parest,allpar,moms),gss,[],[],[],[],lb,ub,[],optcon);
		if i>1
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n \n Mean                    %5.4f  %5.4f %5.4f %5.4f \n St. Dev                  (%5.4f)  (%5.4f) (%5.4f) (%5.4f)\n\n',j,i,mcest(i,:), obj(i),mean(mcest(1:i,:)),var(mcest(1:i,:)).^0.5);
		else
			fprintf('Simulation %d :: Trial %d :: Estimate = (%5.4f,%5.4f,%5.4f,%5.4f) | Obj = %5.4f \n \n ',j,i,mcest(i,:), obj(i));
		end
	end
%  	save(['out1/psimest' int2str(j) '.mat']);
end
%  [parest2 obj2]	=	fminsearch( @(par) psimspecestobj002(par,pcons(1,:),zeros(1,T),zeros(1,T),x,xvar,H,N,b),pargss,opt);

%  [parest1 obj1]

