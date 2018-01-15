%  Copyright 2012 Timothy John Schwuchow
%  psimspeceval002.m		-	Evaluate speculator behavior
clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');

load('outmat/specparam.mat');
j=1;
for j=1:ngrid


	fprintf('Analysis for simulation %d \n __________________________\n',j);
	fprintf('\n Position Volatility Analysis\n\n');
	load(['outmat/pol' int2str(j) '.mat']);
	load(['out1/csvf' int2str(j) '.mat']);
	mpol			=	max(pol(1,:));
	stlist		=	smap(ceil(epts/2),:,:);
	mplist		=	pol(1,:)==mpol;
	bigstates	=	smapinv(mplist,:);
	bigstatemat	=	[umap(j,bigstates(:,2))' dmap(j,bigstates(:,3))'];
	fprintf('\n \n States for largest inventory from zero to %4.3f \n \n',imap(mpol));
	fprintf([' u  |  delta ' repmat('\n  %5.4f  | %5.4f ',1,size(bigstatemat,1)) '\n\n'],bigstatemat');

	polist		=	pol(1,stlist);
	sinvmat		=	[umap(j,smapinv(stlist,2))' dmap(j,smapinv(stlist,3))' imap(polist)'];
%  	fprintf('\n \n States and inventory changes \n \n');
%  	fprintf([' u  |  delta | h' repmat('\n  %5.4f  | %5.4f  || | %5.4f ',1,size(sinvmat,1)) '\n\n'],sinvmat');

%  	polslice		=	pol(1,smap(ceil(epts/2),:,20));
%  	fprintf('Sample inventory jump for eta=28,u=%d / %d, d = 20 :: h %5.4f -> %5.4f \n\n',32,33,imap(polslice(32:33)));
%  	sslice		=	smap(28,32:33,20);
%  	[emap(j,smapinv(sslice,1))' umap(j,smapinv(sslice,2))' dmap(j,smapinv(sslice,3))']
%  	flowprice(:,sslice)

%  	for i=1:2
%  		fprintf('Flow price %d :: %5.4f \n',i,pricespec002(imap(polslice(i+31)),0,emap(j,smapinv(sslice(i),1)),dmap(j,smapinv(sslice(i),3)),umap(j,smapinv(sslice(i),2)),x,b,regrid(j),H,N,bmvargrid(j),cd1grid(j),ce1grid(j)));
%  	end
end

fprintf('\n');
