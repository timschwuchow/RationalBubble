function [V pol] = VFMax002_tx(EV,pricegrid,imap,imapstate,b,tx)
	ipts	=	numel(imap);
	pol	=	zeros(size(EV));
	V		=	zeros(size(EV));
	for j = 1:ipts
		[V(j,:),pol(j,:)]	=	max(( (imap(j) - imapstate) - tx*abs(imap(j) - imapstate)).*pricegrid + b * EV,[],1);
	end