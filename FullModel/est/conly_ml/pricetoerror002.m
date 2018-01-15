%  Copyright 2012 Timothy John Schwuchow
%  pricetoerror002.m 	-	Isolates residual error term (maximum likelihood version)

function [phat d_phat dd_phat]	=	pricetoerror002(p,inv,invexp,b,x,H,N,bmvar,re,cd2)

T			=	numel(p);

p			=	reshape(p,[1 T]);
inv		=	reshape(inv,size(p));
invexp	=	reshape(invexp,size(p));
d_T		=	T-1;

bm			=	norminv(1-(H-inv)/N,0,bmvar);
bmexp		=	norminv(1-(H-invexp)/N,0,bmvar);


phat		=	(p - (x + b*bmexp)/(1-b) - bm)';
d_phat	=	phat(2:T) - cd2*phat(1:T-1);
dd_phat	=	d_phat(2:d_T) - re*d_phat(1:d_T-1);


