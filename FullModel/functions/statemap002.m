%  Copyright 2012 Timothy John Schwuchow
%  statemap002.m		-	Generate state mapping
function [imap emap umap dmap smap smapinv nstate]	=	statemap002(pgrid,vdgrid,cb2grid,cd2grid,bmvargrid,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N)


nstate  		=  epts*upts*dpts;
evec			=	(1:epts);
dvec			=	(1:dpts);
ivec			=	(1:ipts);
uvec			=	(1:upts);
smap			=	reshape((1:nstate),[epts,upts,dpts]);
smapinv		=	[repmat(evec',dpts*upts,1) repmat(reshape(repmat(uvec,epts,1),[upts*epts 1]),dpts,1) reshape(repmat(dvec,epts*upts,1),[nstate 1])];


regrid		=	pgrid(:,3);
vegrid		=	pgrid(:,1);
imap			=	linspace(0,maxshare*H,ipts);
etasd			=	(1/(1-regrid.^2)).^0.5;
emap			=	zeros(numel(regrid),epts);
umap			=	zeros(numel(regrid),upts);
dmap			=	zeros(numel(regrid),dpts);
bmmin			=	zeros(numel(regrid),1);
bmmax			=	bmmin;
dsdmult		=	(1/(1-cd2grid)).^0.5;
for i=1:numel(regrid)
	umap(i,:)	=	linspace(-usd*vegrid(i).^0.5,usd*vegrid(i).^0.5,upts);
	emap(i,:)	=	linspace(-esd*etasd(i)*vegrid(i)^0.5,esd*etasd(i)*vegrid(i)^0.5,epts);
	bmmin(i)		=	norminv(1-H/N,0,bmvargrid(i));
	bmmax(i)		=	norminv(1-(H - maxshare*H)/N,0,bmvargrid(i));
	dmap(i,:)	=	linspace(-dsd*dsdmult(i)*vdgrid(i)^0.5,dsd*dsdmult(i)*vdgrid(i)^0.5 + ((1/(1-cd2grid(i))).^0.5*(bmmax(i) - bmmin(i)))/dsdmult(i),epts);
end
save('states.mat','imap','emap','umap','dmap','smap','smapinv','nstate','-V7.3');
