%  Copyright 2012-2013 Timothy John Schwuchow
%  VFMax002.m		-	Maximize VF
%  Production version

function [V pol] = VFMax002(EV,pricegrid,imap,imapstate,b)
	ipts	=	numel(imap);
	pol	=	zeros(size(EV));
	V		=	zeros(size(EV));
	for j = 1:ipts
		[V(j,:),pol(j,:)]	=	max(( (imap(j) - imapstate)).*pricegrid + b * EV,[],1);
	end