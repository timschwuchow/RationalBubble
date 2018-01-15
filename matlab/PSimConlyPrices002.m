%  Copyright 2012-2013 Elliot Anenberg, Patrick Bayer, James Roberts, and Timothy John Schwuchow
%  PSimConlyPrices002.m		-	Template file (used by SimConlyPrices) - Instances of 'r-e-p-n-g' are replaced by simulation trial number, after which the entire file is copied to a simulation-specific version and run.
%  Current Version

%clear all
%format short g
%path(path,'../../functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load common parameter vectors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('outmat/specparam.mat');

randn('state',seed);
rand('state',seed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for current trial %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcons						=	zeros(nsim,T,8);
for repng=1:8
simstate		=	repng;
x				=	xgrid(repng);
xvar			=	xvargrid(repng);
re			=	regrid(repng);
bmvar		=	bmvargrid(repng);
%varvgrid		=	varvgrid(repng);
cd1			=	cd1grid(repng);
ce1			=	ce1grid(repng);
ve			=	vegrid(repng);
cd2			=	cd2grid(repng);
ce2			=	ce2grid(repng);
cb2			=	cb2grid(repng);
mm          =   mmgrid(repng);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Draw common errors and simulate price sequence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta						=	zeros(nsim,T);
eta							=	zeros(nsim,T);
u							=	normrnd(0,ve^0.5,[nsim T]);
z                           =   normrnd(0,vz^0.5,[nsim T]);

%calculate the marginal bidder given the noise trader error z
 bnew=ones(nsim,T);bold=zeros(nsim,T);iter=0;
 while(max(max(abs(bnew-bold)>.0001)))
     max(abs(bnew-bold))
 bold=bnew;    
 bnew=norminv(1-H-q*(normcdf(bold-z,0,bmvar)-normcdf(bold,0,bmvar)),0,bmvar);   
 iter=iter+1;
 end
 bm=bnew;
 bmexp=ones(nsim,T)*mm;       
        

for t=1:T
	if t > 1
		eta(:,t) 	=	re*eta(:,t-1) + u(:,t);
        
                
      
		delta(:,t)	=	cd2*delta(:,t-1) + ce2*u(:,t-1)+cb2*(bm(:,t-1)-bmexp(:,t-1));				% Delta scripted t=t-1 (t is the delta that determines price in period t)
	else
		eta(:,t) =	u(:,t);
	end
	pcons(:,t,repng)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t),delta(:,t),u(:,t),x,b,re,H,N,bmvar,cd1,ce1,z(:,t),q,mm,bm(:,t),bmexp(:,t));
end

%peff						=	zeros(nsim,T);
%for t=1:T
	%peff(:,t)			=	pricespec002(zeros(nsim,1),zeros(nsim,1),eta(:,t),zeros(nsim,1),u(:,t),x,b,regrid,H,N,bmvargrid,0,0);
%end
%peff is perfectly informed case?  can ignore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save all data and price sequence csv files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save(['outmat/pricesimconly' int2str(repng) '.mat'],'-V7.3');
% 
% datout	=	[repng*ones(T,1) (1:T)' pcons(1,:)' peff(1,:)'];
% csvwrite(['outcsv/PriceSim' int2str(repng) '.csv'],datout);

end
%  quit
