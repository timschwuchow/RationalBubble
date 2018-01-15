%  Copyright 2012 Timothy John Schwuchow
%  nsprobpar002_WV.m		-	Generate speculator transition kernel (parallel)
%  Working version
function [cnsprob]	=	nsprobpar002_WV(ve,re,ce2,cd2,cb2,bmvar,b,H,N,imap,emap,umap,dmap,smapinv,smap,nonaive)

%  Unpack

ipts						=	numel(imap);
epts						=	numel(emap);
upts						=	numel(umap);
dpts						=	numel(dmap);
nstate					=	epts*upts*dpts*ipts;




if nonaive==1
	naiveerr				=	zeros(size(imap));
else
	naiveerr				=	(norminv(1 - (H)/ N,0,bmvar) - norminv(1 - (H - imap(smapinv(:,1))')/N,0,bmvar));
end
t							=	tic;

t							=	tic;
fprintf('Constructing delta mapping\n')
dbase						=	repmat(cd2*dmap(smapinv(:,4))'+ce2*umap(smapinv(:,3))'+cb2*naiveerr,[1 dpts]);
dblock					=	repmat(dmap,[nstate 1]);
[~,ciermap]				=	min(abs(dblock -dbase),[],2);
clear dblock dbase
ciermap					=	reshape(repmat(ciermap',[upts,1]),[nstate*upts 1]);
fprintf('Finished in %5.4f seconds\n',toc(t));
fprintf('Constructing eta mapping\n')
t							=	tic;
eblock					=	repmat(emap,[epts*upts 1]);
[~,enextstate]			=	min(abs(eblock - repmat(re*reshape(repmat(emap,[upts 1]),[epts*upts 1]) + repmat(umap',[epts,1]),[1 epts])),[],2);
enextstate				=	repmat(enextstate,[dpts*upts*ipts 1]);
clear eblock
fprintf('Finished in %5.4f seconds\n',toc(t));
fprintf('Constructing overall mapping\n')
t							=	tic;
unextstate				=	repmat((1:upts)',nstate,1);
inextstate				=	reshape(repmat(smapinv(:,1)',upts,1),[nstate*upts 1]);
spind						=	sub2ind(size(smap),inextstate,enextstate,unextstate,ciermap);
fprintf('Finished in %5.4f seconds\n',toc(t));
fprintf('Constructing final mapping\n')
t							=	tic;
clear inextstate enextstate ciermap unextstate smap smapinv

punorm					=	normpdf(umap',0,ve^0.5);
punorm					=	repmat(punorm/sum(punorm),[nstate 1]);
sind						=	reshape(repmat((1:nstate),upts,1),[nstate*upts 1]);
fprintf('Finished in %5.4f seconds\n',toc(t));
fprintf('Constructing transition kernel\n')
t							=	tic;
cnsprob					=	sparse(spind,sind,punorm,nstate,nstate);
fprintf('Finished in %5.4f seconds\n',toc(t));



