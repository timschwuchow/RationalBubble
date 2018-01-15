%  Copyright 2012 Timothy John Schwuchow
%  nsprobpar002.m		-	Generate speculator transition kernel (parallel)
%  Production version
function [cnsprobpart]	=	nsprobpar002(ve,re,ce2,cd2,cb2,bmvar,b,H,N,imap,emap,umap,dmap,smapinv,smap,SimNum,nonaive)

%  Unpack

ipts						=	numel(imap);
epts						=	numel(emap);
upts						=	numel(umap);
dpts						=	numel(dmap);
nstate					=	epts*upts*dpts;

dblock					=	repmat(dmap,[nstate 1]);

if nonaive==1
	naiveerr				=	zeros(size(imap));
else
	naiveerr				=	(norminv(1 - H/ N,0,bmvar) - norminv(1 - (H - imap(SimNum))/N,0,bmvar));
end


dbase						=	repmat(cd2*dmap(smapinv(:,3))'+ce2*umap(smapinv(:,2))'+cb2*naiveerr,[1 dpts]);
[~,ciermap]				=	min(abs(dblock -dbase),[],2);
clear dbase dblock naiveerr
ciermap					=	reshape(repmat(ciermap',[upts,1]),[nstate*upts 1]);
eblock					=	repmat(emap,[epts*upts 1]);
[~,enextstate]			=	min(abs(eblock - repmat(re*reshape(repmat(emap,[upts 1]),[epts*upts 1]) + repmat(umap',[epts,1]),[1 epts])),[],2);
clear eblock
enextstate				=	repmat(enextstate,[dpts*upts 1]);
unextstate				=	repmat((1:upts)',nstate,1);
spind						=	sub2ind(size(smap),enextstate,unextstate,ciermap);
clear enextstate ciermap unextstate  smap smapinv
punorm					=	normpdf(umap',0,ve^0.5);
punorm					=	repmat(punorm/sum(punorm),[nstate 1]);
sind						=	reshape(repmat((1:nstate),upts,1),[nstate*upts 1]);
cnsprobpart				=	sparse(spind,sind,punorm,nstate,nstate);



