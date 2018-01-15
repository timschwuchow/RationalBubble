imapstate	=	repmat(imap',1,nstate);
emapstate	=	repmat(emap(48,smapinv(:,1)),[ipts 1]);
umapstate	=	repmat(umap(48,smapinv(:,2)),[ipts 1]);
dmapstate	=	repmat(dmap(48,smapinv(:,3)),[ipts 1]);
pricegrid	=	vfprice002(imap,zeros(size(imap)),emapstate(1,:),dmapstate(1,:),umapstate(1,:),xgrid(48),b,regrid(48),H,N,bmvargrid(48),cd1grid(48),ce1grid(48));
stateextend	=	reshape(repmat((1:nstate),[ipts 1]),[nstate*ipts 1]);
for j=1:ipts
	t1						=	tic;
	[cnsprob{j}]		=	nsprobpar002(pgrid(48,1),pgrid(48,3),ce2grid(48),cd2grid(48),cb2grid(48),bmvargrid(48),b,H,N,imap,emap(48,:),umap(48,:),dmap(48,:),smapinv,smap,j,nonaive);
	fprintf('Transition matrix %d solved in %5.4f minutes. \n',j,toc(t1)/60);
end
price	=	zeros(ipts,2);
flow	=	price;
ev		=	price;
for i=1:ipts
	for n=1:2
		price(i,n)	=	pricegrid(i,s(n));
		flow(i,n)	=	(imap(1) - imap(i))*pricegrid(i,s(n));
		ev(i,n)		=	b*V(i,:)*cnsprob{i}(:,s(n));
	end
end
[price imap']
FlowFix = zeros(size(flow));
for i=1:ipts; for n=1:2; FlowFix(i,n) =	flow(i,n) - norminv(1-(H-imap(i))/N,0,bmvargrid(48));end;end;
CSVF			=	flow + ev
fig=figure;h1=line(imap(10:30),ev(10:30,1),'Color','r');hold all;h2=line(imap(10:30),-flow(10:30,1),'Color','b','Parent',gca);print(fig,'-djpeg','EVFlow1,jpg');close(fig);fig = figure;h3=line(imap(10:30),ev(10:30,2),'Color','r');hold all;h4=line(imap(10:30),-flow(10:30,2),'Color','b','Parent',gca);print(fig,'-djpeg','EVFlow2.jpg');close(fig);