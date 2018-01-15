%  Copyright 2012 Timothy John Schwuchow
%  psimcsvf002.m		-	Create csvfs

format short g
load('outmat/specparam.mat');
path(path,'../../functions');
path(path,'../../functions/tests');

imapstate	=	repmat(imap',1,nstate);

for n=1:ngrid
	fprintf('Graphing trial %d \n \n', n);
	load(['outmat/pol' int2str(n) '.mat']);
	cellnsprob	=	cell(ipts,1);
	fprintf('Loading transition kernel\n');
	for j=1:ipts
		load(['outmat/trans' int2str(j) '_' int2str(n) '.mat']);
		cellnsprob{j}		=	cnsprob;
	end
	cnsprob	=	cellnsprob;
	clear cellnsprob
	EV					=	zeros(ipts,nstate);
	fprintf('Computing EV \n');
	for j = 1:ipts
		EV(j,:)		=	V(j,:)*cnsprob{j};
	end
	emapstate	=	repmat(emap(n,smapinv(:,1)),[ipts 1]);
	dmapstate	=	repmat(dmap(n,smapinv(:,3)),[ipts 1]);
	umapstate	=	repmat(umap(n,smapinv(:,2)),[ipts 1]);
	dslfd
	csvf			=	cell(numel([1:5:ipts]),1);
	fprintf('Computing CSVF\n \n');
	flowprice	=	vfprice002(imap,zeros(size(imap)),emapstate(1,:),dmapstate(1,:),umapstate(1,:),x,b,regrid(n),H,N,bmvargrid(n),cd1grid(n),ce1grid(n),nstate);
	k				=	1;
	for j = 1:5:ipts
		csvf{k}	=	(imap(j) - imapstate).*flowprice + b * EV;
		k			=	k+1;
	end
	save(['out1/csvf' int2str(n) '.mat'],'flowprice','csvf','EV','-v7.3');
	clear csvf EV cnsprob
end


quit