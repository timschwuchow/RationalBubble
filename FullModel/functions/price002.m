% Copyright 2012 Patrick Bayer, James Roberts, and Timothy Schwuchow
% price002.m	-	Compute prices

function p	=	price002(inv,invexp,eta,delta,inno,x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid)

bm				=	zeros(size(regrid));
nsim			=	numel(regrid);
for i=1:nsim
	bm		=	norminv(1 - (H-inv(i))/ N,0,bmvargrid(i));
	bmexp		=	norminv(1 - (H-invexp(i))/ N,0,bmvargrid(i));
end
p		=	(x + b * bmexp) / (1 - b) + bm + (eta + b*regrid.*(cd1grid.*delta + ce1grid.*inno ))./ (1 - b * regrid);