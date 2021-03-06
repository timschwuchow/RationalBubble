function V = CSVFHowardImp002_WV(pol,stateextend,V,cnsprob,pricegrid,b,imap,imapstate,maxiter,tolcrit)
path(path,'../../functions');

pol		=	cell2mat(pol);
V			=	cell2mat(V);

ipts		=	size(pol,1);
nstate	=	size(pol,2);
polshape	=	reshape(pol,[numel(pol),1]);
polind	=	reshape(sub2ind(size(V),polshape,stateextend),size(pol));
polprice	=	reshape(pricegrid(sub2ind(size(pricegrid),polshape,stateextend)),size(V));
polimap	=	imap(pol);
tol		=	1;
iter		=	0;
Vo			=	V;
while ( (tol > tolcrit) & (iter < maxiter))
	iter	=	iter+1;
	EV		=	EVEval002(V,cnsprob);
	V		=	(imapstate - polimap).*polprice + b*EV(polind);
	tol	=	max(max(abs(V-Vo)));
	Vo		=	V;
end
%  fprintf('Howard iterations %d \n',iter);
V			=	mat2cell(V,ones(ipts,1),nstate);