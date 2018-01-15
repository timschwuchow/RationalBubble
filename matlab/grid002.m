%  Copyright 2012-2013 Timothy John Schwuchow
%  grid002.m		-	Generate parameter grid
%  Production Version
function [vegrid vggrid regrid rggrid xgrid xvargrid pgrid nsim]	=	grid002(ve,vg,re,rg,x,xvar,b,vdiag,q,vz);

if vdiag~=1
	rggrid	=	repmat(rg',numel(vg)*numel(ve)*numel(re)*numel(x)*numel(xvar),1);
	regrid	=	repmat(reshape(repmat(re,numel(rg),1),[numel(rg)*numel(re) 1]),numel(vg)*numel(ve)*numel(x)*numel(xvar),1);
	vggrid	=	repmat(reshape(repmat(vg,numel(re)*numel(rg),1),[numel(rg)*numel(re)*numel(vg) 1]),numel(ve)*numel(x)*numel(xvar),1);
	vegrid	=	repmat(reshape(repmat(ve,numel(re)*numel(rg)*numel(vg),1),[numel(rg)*numel(re)*numel(vg)*numel(ve) 1]),numel(x)*numel(xvar),1);
	xgrid		=	repmat(reshape(repmat(x,numel(re)*numel(rg)*numel(vg)*numel(ve),1),[numel(rg)*numel(re)*numel(vg)*numel(ve)*numel(x) 1]),numel(xvar),1);
	xvargrid	=	reshape(repmat(xvar,numel(re)*numel(rg)*numel(vg)*numel(ve)*numel(x),1),[numel(rg)*numel(re)*numel(vg)*numel(ve)*numel(x)*numel(xvar) 1]);
	nsim		=	numel(vegrid);
	pgrid		=	[vegrid vggrid regrid rggrid xgrid xvargrid];
else
	minv 	= 	min([numel(ve),numel(vg),numel(re),numel(rg),numel(x),numel(xvar)]);
	rggrid	=	rg(1:minv)';
	regrid	=	re(1:minv)';
	vggrid	=	vg(1:minv)';
	vegrid	=	ve(1:minv)';
	xgrid		=	x(1:minv)';
	xvargrid	=	xvar(1:minv)';
	nsim		=	numel(vegrid);
	pgrid		=	[vegrid vggrid regrid rggrid xgrid xvargrid];
end