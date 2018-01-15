%  Copyright 2012 Timothy John Schwuchow
%  PSimTrans002.m		-	Generate transition kernel for speculator transition (parallel)
%  Production version
clear all
format short g
path(path,'../../functions');

load('outmat/specparam.mat');


cnsprob				=	cell(ipts,1);
for j=1:ipts
	[cnsprob{j}]		=	nsprobpar002(pgrid(repng,1),pgrid(repng,3),ce2grid(repng),cd2grid(repng),cb2grid(repng),bmvargrid(repng),b,H,N,imap,emap(repng,:),umap(repng,:),dmap(repng,:),smapinv,smap,j,nonaive);
end
save(['outmat/trans' int2str(repng) '.mat'],'cnsprob','-V7.3');

quit

