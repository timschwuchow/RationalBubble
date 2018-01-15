%  Copyright 2012-2013 Timothy John Schwuchow
%  CSVFEVEval002.m		-	Evaluate Expected VF for CSVF formulation
%  Production version

function EV = CSVFEVEval002(V,cnsprob)
	ipts				=	size(V,1);
	EV					=	zeros(ipts,numel(V{1}));
	for k=1:ipts
		EV(k,:)		=	V{k}*cnsprob{k};
	end
end