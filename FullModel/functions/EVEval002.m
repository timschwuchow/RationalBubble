%  Copyright 2012-2013 Timothy John Schwuchow
%  EVEval002.m		-	Evaluate Expected VF
%  Production version

function EV = EVEval002(V,cnsprob)
	ipts				=	size(V,1);
	EV					=	zeros(size(V));
	for k=1:ipts
		EV(k,:)		=	V(k,:)*cnsprob{k};
	end
end