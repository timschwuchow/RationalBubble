function V = HowardImp002(pol,stateextend,V,cnsprob,pricegrid,b,imap,imapstate,howimp)

ipts		=	size(pol,1);
nstate	=	size(pol,2);
polshape	=	reshape(pol,[numel(pol),1]);
polind	=	reshape(sub2ind(size(V),polshape,stateextend),size(pol));
for j= 1:howimp
	EV		=	EVEval002(V,cnsprob);
	V		=	(imapstate - imap(pol)).*reshape(pricegrid(sub2ind(size(pricegrid),polshape,stateextend)),size(V)) + b*EV(polind);
end