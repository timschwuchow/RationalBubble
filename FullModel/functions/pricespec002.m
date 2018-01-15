% Copyright 2012 Patrick Bayer, James Roberts, and Timothy Schwuchow
% pricespec002.m	-	Compute prices (uses change of notation to correspond to brstheory_eq001.pdf theory results

function p	=	pricespec002(inv,invexp,eta,delta,inno,x,b,re,H,N,bmvar,cd1,ce1)


bm			=	norminv(1 - (H-inv)/ N,0,bmvar);
bmexp		=	norminv(1 - (H-invexp)/ N,0,bmvar);

p		=	(x + b * bmexp) / (1 - b) + bm + eta ./ (1 - b * re) +  cd1.*delta + ce1.*inno;