%  Copyright 2012 Timothy John Schwuchow
%  vfsol002.m		-	Solve speculator value/policy function (parallel)
%  Working version

function [pol V]	=	vfsolpar002_WV(re,cd1,ce1,bmvar,b,H,N,x,tolcrit,maxiter,cnsprob,imap,emap,dmap,umap,smap,smapinv,howimp,howlag)


tstart		=	tic;
nstate		=	numel(emap)*numel(umap)*numel(dmap)*numel(imap);
ipts			=	numel(imap);

imap(smapinv(:,1))'*vfprice002(zeros(nstate,1),zeros(nstate,1),emap(smapinv(:,2))',dmap(smapinv(:,4))',umap(smapinv(:,3))',x,b,re,H,N,bmvar,cd1,ce1);
%  V					=	repmat(imap',[1  nstate]).*vfprice002(zeros(size(imap)),zeros(size(imap)),emap(smapinv(:,1)),dmap(smapinv(:,3)),umap(smapinv(:,2)),x,b,re,H,N,bmvar,cd1,ce1);
pol				=	smapinv(:,1);
%  pol				=	repmat((1:ipts)',1,nstate);
tol				=	1;
iter				=	0;
Vo					=	V;
EV					=	V;
%  imapstate	=	repmat(imap',1,nstate);
%  emapstate	=	repmat(emap(smapinv(:,1)),[ipts 1]);
%  umapstate	=	repmat(umap(smapinv(:,2)),[ipts 1]);
%  dmapstate	=	repmat(dmap(smapinv(:,3)),[ipts 1]);
%  pricegrid	=	vfprice002(imap,zeros(size(imap)),emapstate(1,:),dmapstate(1,:),umapstate(1,:),x,b,re,H,N,bmvar,cd1,ce1);
%  stateextend	=	reshape(repmat((1:nstate),[ipts 1]),[nstate*ipts 1]);

while ( ( tol > tolcrit) & (iter < maxiter) )
	t1			=	tic;
	iter		=	iter + 1;
	if iter>howlag
		polshape	=	reshape(pol,[nstate*ipts,1]);
		polind	=	reshape(sub2ind(size(EV),polshape,stateextend),[ipts,nstate]);
		for k=1:ipts
			EV(k,:)		=	Vo(k,:)*cnsprob{k};
		end
		for j= 1:howimp
			V		=	(imapstate - imap(pol)).*reshape(pricegrid(sub2ind(size(pricegrid),polshape,stateextend)),size(V)) + b*EV(polind);
			for k=1:ipts
				EV(k,:)		=	V(k,:)*cnsprob{k};
			end
		end
		for j = 1:ipts
			[V(j,:),pol(j,:)]	=	max((imap(j) - imapstate).*pricegrid + b * EV,[],1);
		end
	else
		for j = 1:ipts
			EV(j,:)		=	Vo(j,:)*cnsprob{j};
		end
		for j = 1:ipts
			[V(j,:),pol(j,:)]	=	max((imap(j) - imapstate).*pricegrid + b * EV,[],1);
		end
	end
	tol			=	max(max(abs((V - Vo))));
	fprintf('iter %d -- tol %5.4g -- in %5.4f minutes --  (%5.4f hours cumulatively) \n', iter, tol, toc(t1)/60,toc(tstart)/60/60);
	Vo			=	V;
end

fprintf('Policy function solved in %d iterations after %5.3f minutes\n',iter,toc(tstart)/60);

