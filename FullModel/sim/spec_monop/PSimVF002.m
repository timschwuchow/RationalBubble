%  Copyright 2012 Timothy John Schwuchow
%  PSimVF002.m		-	Solve VF and simulate (parallel)
%  Production version
clear all
format short g
path(path,'../../functions');

load('outmat/specparam.mat');
cnsprob				=	cell(ipts,1);
t0							=	tic;
for j=1:ipts
	t1						=	tic;
	[cnsprob{j}]		=	nsprobpar002(pgrid(repng,1),pgrid(repng,3),ce2grid(repng),cd2grid(repng),cb2grid(repng),bmvargrid(repng),b,H,N,imap,emap(repng,:),umap(repng,:),dmap(repng,:),smapinv,smap,j,nonaive);
	fprintf('Transition matrix %d solved in %5.4f minutes. \n',j,toc(t1)/60);
end
fprintf('\n ************************ \nTransition kernel solved in %5.4f seconds \n\n',toc(t0));
t0				=	tic;
V				=	repmat(imap',[1  nstate]).*vfprice002(zeros(size(imap)),zeros(size(imap)),emap(repng,smapinv(:,1)),dmap(repng,smapinv(:,3)),umap(repng,smapinv(:,2)),xgrid(repng),b,regrid(repng),H,N,bmvargrid(repng),cd1grid(repng),ce1grid(repng));
[pol V]		=	vfsolpar002(regrid(repng),cd1grid(repng),ce1grid(repng),bmvargrid(repng),b,H,N,xgrid(repng),vftolcrit,vfmaxiter,cnsprob,imap,emap(repng,:),dmap(repng,:),umap(repng,:),smap,smapinv,howimp,howlag,V);
clear cnsprob
save(['outmat/pol' int2str(repng) '.mat'],'pol','V','-V7.3');


quit

