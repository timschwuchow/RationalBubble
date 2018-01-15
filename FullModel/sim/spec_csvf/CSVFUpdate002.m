%  Copyright 2012-2013 Timothy John Schwuchow
%  CSVFUpdate002.m		-	Solve for CSVF,VF,Pol
%  Production version
function [CSVFCurrent VCurrent PolCurrent]	=	 CSVFUpdate002(CSVFLast,VLast,PolLast,ve,re,cd1,ce1,cd2,ce2,cb2,bmvar,b,H,N,x,tolcrit,maxiter,imap,emap,dmap,umap,smap,smapinv,howimp,howlag,nonaive);
path(path,'../../functions');
nstate		=	numel(emap)*numel(umap)*numel(dmap);
ipts			=	numel(imap);
cnsprob		=	cell(ipts,1);
t1						=	tic;
for j=1:ipts
	[cnsprob{j}]		=	nsprobpar002(ve,re,ce2,cd2,cb2,bmvar,b,H,N,imap,emap,umap,dmap,smapinv,smap,j,nonaive);
	if mod(j,5)==0
		fprintf('Transition matrix %d of %d solved in %5.4f seconds. \n',j,ipts,toc(t1));
		t1						=	tic;
	end
end



imapstate	=	repmat(imap',1,nstate);
emapstate	=	repmat(emap(smapinv(:,1)),[ipts 1]);
umapstate	=	repmat(umap(smapinv(:,2)),[ipts 1]);
dmapstate	=	repmat(dmap(smapinv(:,3)),[ipts 1]);
pricegrid	=	vfprice002(imap,zeros(size(imap)),emapstate(1,:),dmapstate(1,:),umapstate(1,:),x,b,re,H,N,bmvar,cd1,ce1);
stateextend	=	reshape(repmat((1:nstate),[ipts 1]),[nstate*ipts 1]);

tol			=	1;
iter			=	0;



while ( ( tol > tolcrit) & (iter < maxiter) )
	t1												=	tic;
	iter											=	iter + 1;
	EV												=	CSVFEVEval002(VLast,cnsprob);
	[~, VCurrent, PolCurrent]	=	CSVFMax002(EV,pricegrid,imap,imapstate,b);
	if iter>=howlag
		VCurrent		=	CSVFHowardImp002_WV(PolCurrent,stateextend,VCurrent,cnsprob,pricegrid,b,imap,imapstate,howimp,1e-5);
	end

	tol	=	max(max(abs(cell2mat(VCurrent) - cell2mat(VLast))));

	fprintf('iter %d -- tol %5.4g -- in %5.4f seconds \n', iter, tol, toc(t1));
	VLast			=	VCurrent;
end
[CSVFCurrent VCurrent PolCurrent]	=	CSVFMax002(EV,pricegrid,imap,imapstate,b);
