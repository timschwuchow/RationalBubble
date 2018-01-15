% Copyright 2012 Patrick Bayer, James Roberts, and Timothy Schwuchow
% pricespec002.m	-	Compute prices (uses change of notation to correspond to brstheory_eq001.pdf theory results

function p	=	pricespec002(inv,invexp,eta,delta,inno,x,b,re,H,N,bmvar,cd1,ce1,z,q,mm,bm,bmexp)

% bnew=1;bold=0;iter=0;
% while(abs(bnew-bold)>.001)
%     bold=bnew;    
% bnew=norminv(1-H-q*(normcdf(bold-z,0,bmvar)-normcdf(bold,0,bmvar)),0,bmvar);   
% iter=iter+1;
% end
% 
% 
% bm=bnew;
% %bm			=	norminv(1 - (H-inv)/ N,0,bmvar);
% %bmexp		=	norminv(1 - (H-invexp)/ N,0,bmvar);
% bmexp=mm;

p		=	(x + b .* bmexp) / (1 - b) + bm + eta ./ (1 - b * re) +  cd1.*delta + ce1.*inno;

%inv is the speculator inventory;
%bmexp is the expected speculator inventory