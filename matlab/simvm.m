function [mm,vm1]		=	simvm(bmvargrid,H,vz,q)

%bmnvargrid=bmvargrid;
%for i=1:length(bmvargrid)
randn('state',8);
T=5000;
zeta=randn(T,1)*sqrt(vz);
bfinal=zeros(T,1);
%for t=1:T;
bnew=ones(T,1);bold=zeros(T,1);iter=0;
while(abs(bnew-bold)>.0001)
    bold=bnew; 
    bnew=norminv(1-H-q*(normcdf(bold-zeta,0,bmvargrid)-normcdf(bold,0,bmvargrid)),0,bmvargrid);    
    iter=iter+1;
end

bfinal=bnew;
%end

[mean(bfinal),std(bfinal)];
vm1=var(bfinal);mm=mean(bfinal);
%end

%[f,xi]=ksdensity(bfinal);
%plot(xi,f);