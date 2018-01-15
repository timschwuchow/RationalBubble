function [CSVFCurrent VCurrent PolCurrent]	=	CSVFIter002(VLast,cnsprob,pricegrid,imap,imapstate,b,flow);

	for i=1:ipts
		flow{i} + repmat(VLast{i},ipts,1)*cnsprob{i}
	end
	EV												=	CSVFEVEval002(VLast,cnsprob);
