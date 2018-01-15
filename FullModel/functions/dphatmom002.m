%  Copyright 2012 Timothy John Schwuchow
%  dphatmom002.m	-	Generates theoretical moments of dphat

function mom	=	dphatmom002(momnum,pgrid,vdgrid,ce1grid,ce2grid,cd1grid,cd2grid,b)

vegrid	=	pgrid(:,1);
regrid	=	pgrid(:,3);
nsim		=	numel(regrid);
mom		=	zeros(nsim,1);
cut		=	(1+b*regrid.*ce1grid);
cud		=	(b.*regrid.*cd1grid.*(cd2grid - regrid));
cutm1		=	(b*regrid.*(cd1grid.*ce2grid - regrid.*ce1grid));
if momnum==0
	mom 	=	cut.^2.*vegrid + cud.^2.*vdgrid + cutm1.^2.*vegrid;
elseif momnum==1
	mom 	=	cut.*cutm1.*vegrid + cud.*(cud.*cd2grid.*vdgrid + ce2grid.*cutm1.*vegrid);
else
	mom	=	cd2grid.^(momnum-2).*ce2grid.*cud.*cut.*vegrid + cd2grid.^(momnum-1).*cud.*ce2grid.*cutm1.*vegrid + cd2grid.^momnum.*cud.^2.*vdgrid;
end