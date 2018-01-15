% Copyright 2012 Patrick Bayer, James Roberts, and Timothy Schwuchow
% vfprice002.m	-	Compute prices in value function solution (for all states)
% Production version
function p	=	vfprice002(inv,invexp,eta,delta,inno,x,b,re,H,N,bmvar,cd1,ce1)


if sum(size(eta)~=size(inv)) > 0
	ninv		=	numel(inv);
	neta		=	numel(eta);
	inv		=	repmat(reshape(inv,[ninv,1]),[1 neta]);
	invexp	=	repmat(reshape(invexp,[ninv 1]),[1 neta]);
	eta		=	repmat(reshape(eta,[1 neta]),[ninv 1]);
	delta		=	repmat(reshape(delta,[1 neta]),[ninv 1]);
	inno		=	repmat(reshape(inno,[1 neta]),[ninv 1]);
end

bm				=	norminv(1 - (H-inv)/ N,0,bmvar);
bmexp			=	norminv(1 - (H-invexp)/ N,0,bmvar);

p				=	(x + b * bmexp) / (1 - b) + bm + eta/(1-b*re) + cd1*delta + ce1*inno;

