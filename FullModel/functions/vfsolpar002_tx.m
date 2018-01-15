%  Copyright 2012 Timothy John Schwuchow
%  vfsol002_tx.m		-	Solve speculator value/policy function (parallel) with transaction cost
%  Production version

function [pol V]	=	vfsolpar002(re,cd1,ce1,bmvar,b,H,N,x,tolcrit,maxiter,cnsprob,imap,emap,dmap,umap,smap,smapinv,howimp,howlag,Vinit,tx)


tstart		=	tic;
nstate		=	numel(emap)*numel(umap)*numel(dmap);
ipts			=	numel(imap);

V					=	Vinit;
pol				=	repmat((1:ipts)',1,nstate);
tol				=	1;
iter				=	0;
Vo					=	V;
EV					=	V;
imapstate	=	repmat(imap',1,nstate);
emapstate	=	repmat(emap(smapinv(:,1)),[ipts 1]);
umapstate	=	repmat(umap(smapinv(:,2)),[ipts 1]);
dmapstate	=	repmat(dmap(smapinv(:,3)),[ipts 1]);
pricegrid	=	vfprice002(imap,zeros(size(imap)),emapstate(1,:),dmapstate(1,:),umapstate(1,:),x,b,re,H,N,bmvar,cd1,ce1);
stateextend	=	reshape(repmat((1:nstate),[ipts 1]),[nstate*ipts 1]);

while ( ( tol > tolcrit) & (iter < maxiter) )
	t1			=	tic;
	iter		=	iter + 1;
	EV			=	EVEval002(Vo,cnsprob);
	[V pol]	=	VFMax002_tx(EV,pricegrid,imap,imapstate,b,tx);
	if iter>=howlag
		V		=	HowardImp002_tx(pol,stateextend,V,cnsprob,pricegrid,b,imap,imapstate,howimp,tx);
	end
	tol			=	max(max(abs((V - Vo))));
	fprintf('iter %d -- tol %5.4g -- in %5.4f minutes --  (%5.4f hours cumulatively) \n', iter, tol, toc(t1)/60,toc(tstart)/60/60);
	Vo			=	V;
end

fprintf('Policy function solved in %d iterations after %5.3f hours\n',iter,toc(tstart)/60/60);

