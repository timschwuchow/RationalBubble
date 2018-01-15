%  Copyright 2012 Timothy John Schwuchow
%  innogen002.m	-	Generate common shock innovations

function inno 	=	innogen002(pgrid,T)

vegrid	=	pgrid(:,1);
nsim		=	numel(vegrid);
inno		=	zeros(nsim,T);

for i=1:nsim
	inno(i,:)	=	normrnd(0,vegrid(i)^0.5,[1 T]);
end
