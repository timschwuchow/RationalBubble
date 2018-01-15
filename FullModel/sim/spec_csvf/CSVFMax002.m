%  Copyright 2012-2013 Timothy John Schwuchow
%  CSVFMax002.m		-	Maximize VF
%  Production version

function [CSVF V pol] = CSVFMax002(EV,pricegrid,imap,imapstate,b)
	ipts		=	numel(imap);
	nstate	=	size(EV,2);
	pol		=	cell(ipts,1);
	V			=	cell(ipts,1);
	CSVF		=	cell(ipts,1);
	for i = 1:ipts
		CSVF{i}			=	(imap(i) - imapstate).*pricegrid + b*EV;
		[V{i},pol{i}]	=	max(CSVF{i},[],1);
	end

