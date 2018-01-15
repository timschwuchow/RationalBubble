%  Copyright 2012 Timothy John Schwuchow
%  vfsol002.m		-	Solve speculator value/policy function

function [pol2]	=	vfsol002(pg,vd,ce2,cd2,cb2,bmvar,b,H,N,x,H,N,tolcrit,maxiter)


re							=	pg(3);

%  [imap emap umap dmap smap smapinv nstate]	=	statemap002(pg,vd,cb2,cd2,bmvar,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);
%  [cnsprob,ciermap]										=	nsprob002(pg,vd,ce2,cd2,cb2,bmvar,b,H,N,imap,emap,ermap)
load('states.mat');
load('trans.mat');


V2					=	repmat(imap',[1 nstate]).*vfprice002(zeros(size(imap)),zeros(size(imap)),emap(smapinv(:,1)),dmap(smapinv(:,3)),umap(smapinv(:,2)),x,b,re,H,N,bmvar);

pol2				=	repmat(imap',1,nstate);
tol2				=	1;
iter2				=	0;
Vo				=	V2;
while ( ( tol2 > tolcrit) & (iter2 < maxiter) )
	iter2		=	iter2 + 1;
	EV			=	zeros(ipts,nstate);
	for j = 1:ipts
		EV(j,:)		=	Vo(j,:)*cnsprob{j};
	end
	for j = 1:ipts
		[V2(j,:),pol2(j,:)]	=	max(repmat(imap(j) - imap,1,nstate).*v2price002(imap,zeros(size(imap)),emap(smapinv(:,1)),dmap(smapinv(:,3)),umap(smapinv(:,2)),x,b,re,H,N,bmvar)' + b * EV);
	end
	if iter2>2
		pol2ind	=	reshape(sub2ind(size(EV),reshape(pol2,[nstate*ipts,1]),reshape(repmat((1:nstate),ipts,1),[nstate*ipts,1])),[ipts,nstate]);
		for j= 1:howimp
			for k=1:ipts
				EV(k,:)		=	V2(k,:)*cnsprob{k};
			end
			V2		=	(repmat(imap', [1 nstate]) - imap(pol2)).*vfprice002(imap(pol2),zeros(size(imap)),repmat(emap(smapinv(:,1)),ipts,1),repmat(dmap(smapinv(:,3)),ipts,1),repmat(umap(smapinv(:,2)),ipts,1),x,b,re,H,N,bmvar) + b * EV(pol2ind);
		end
		for j = 1:ipts
			EV(j,:)		=	V2(j,:)*cnsprob{j};
		end
		for j = 1:ipts
			[V2(j,:),pol2(j,:)]	=	max(repmat(imap(j) - imap,1,nstate).*priceprod(repmat(imap',nstate,1),zeros(nstate,ipts),emap(smapinv(:,1)),ermap(smapinv(:,2)),x,b,re,H,N,bmvar)' + b * EV);
		end
	end
	tol2			=	max(max(abs((V2 - Vo))));
	fprintf('iter %d tol %6.4g \n', iter2, tol2)
	Vo			=	V2;
end

