%  Copyright 2012 Timothy John Schwuchow
%  nsprob002.m		-	Generate speculator transition kernel

function [cnsprob,ciermap]	=	nsprob002(pg,vd,ce2,cd2,cb2,bmvar,b,H,N,ipts,epts,upts,dpts,j)

%  Unpack
ve			=	pg(1);
vg			=	pg(2);
re			=	pg(3);
rg			=	pg(4);



%% State mapping %%

load('states.mat');

eblock		=	repmat(emap,[nstate*upts 1]);
dblock		=	repmat(dmap,[nstate*upts 1]);
naiveerr		=	(norminv(1 - (H-imap)/ N,0,bmvar) - norminv(1 - H/ N,0,bmvar));
sind			=	reshape(repmat((1:nstate),upts,1),[nstate*upts 1]);
dbase			=	repmat(cd2*dmap(smapinv(sind,3))' + ce2*umap(smapinv(sind,2))',[1 dpts]);
unextstate	=	repmat((1:upts)',nstate,1);
probu			=	normpdf(umap(unextstate)',0,ve^0.5);
ciermap		=	cell(ipts,1);
cnsprob		=	cell(ipts,1);
clear uvec
fprintf('ipt = %d \n',j);
[~,ciermap{j}]			=	min(abs(dblock - (dbase + cb2*naiveerr(j))),[],2);
[~,enextstate]			=	min(abs(eblock - repmat(re*emap(smapinv(sind,1))' + umap(unextstate)',[1 epts])),[],2);
spind						=	sub2ind(size(smap),enextstate,unextstate,ciermap{j});
clear enextstate
cnsprob{j}				=	sparse(spind,sind,probu,nstate,nstate);
clear spind
cnsprob{j}				=	cnsprob{j}./ repmat(sum(cnsprob{j},1),[nstate,1]);

save(['trans' int2str(j) '.mat'],'ciermap','cnsprob','-V7.3');
