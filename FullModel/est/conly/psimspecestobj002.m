%  Copyright 2012 Timothy John Schwuchow
%  psimspecest002.m		- Objective function for estimated speculation model

function obj	= 	psimspecestobj002(par,parset,price,inv,invexp,x,xvar,H,N,b,parest,allpar,moms)


parvec			=	zeros(1,numel(allpar));
ct					=	1;
for i=1:numel(allpar)
	k	=	strmatch(allpar{i},parest);
	if isempty(k)
		parvec(i)	=	parset(i);
	else
		parvec(i)	=	par(ct);
		ct				=	ct+1;
	end
end
%  Unpack
ve		=	parvec(1);
vg		=	parvec(2);
re		=	parvec(3);
rg		=	parvec(4);


%  Generate constants
[vd vi ce1 cg1 cd1 ci1 ce2 cg2 cd2 ci2 cb2 bmvar varv]		=	bayes002_t1(parvec,b,xvar);

[phat dprice ddprice]		=	pricetoerror002(price,inv,invexp,b,x,H,N,bmvar,re,cd2);
T				=	numel(ddprice);
errcfs		=	[ (1/(1-b*re) + b*re/(1-b*re)*ce1) , -cd2 - b*re*cd2*ce1 + b*re*cd1*ce2 - b*re^2*ce1, b*re^2*(ce1*cd2 - cd1*ce2)];


errcfs		=	[errcfs, zeros(1,moms)];

tmoms			=	[0];
for i=0:moms
	tmoms		=	[tmoms sum(ve*errcfs(1:3).*errcfs(1+i:3+i))];
end



smoms			=	[mean(ddprice)];

for i=0:moms
	cf			=	cov(ddprice(1:T-i),ddprice(1+i:T));
	smoms		=	[smoms, cf(2)];
end

obj			=	mean((tmoms-smoms).^2);

