%  Copyright 2012 Timothy John Schwuchow
%  ErrorVector002.m 	-	Solves for errors in terms of the first two errors and model parameters

function UVec = ErrorVector002(ddprice,u1,u2,b,re,ce1,ce2,cd1,cd2)
psi0 			=	ce1 + 1/(1-b*re);
psi1			=	cd1*ce2 - ce1*cd2 - cd2/(1-b*re) - re*ce1;
psi2			=	re*(ce1*cd2 - cd1*ce2);


ddprice(1) 	= 	ddprice(1) - psi1*u2 - psi2*u1;
ddprice(2)	=	ddprice(2) - psi2*u2;

T				=	numel(ddprice);

Psi													=	eye(T)*psi0;
Psi(sub2ind(size(Psi),[2:T],[1:T-1]))		=	psi1;
Psi(sub2ind(size(Psi),[3:T],[1:T-2]))		=	psi2;

%  Psi			=	sparse(Psi);
%  opts.LT		=	true;

%  UVec			=	linsolve(Psi,ddprice',opts);
%  UVec				=	inv(Psi)*ddprice';
UVec			=	Psi^(-1)*ddprice;

