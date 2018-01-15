

clear all
ngrid					=	50;
H						=	0.3;
N						=	1;
maxshare				=	0.999;
[bmgrid invgrid]	=	meshgrid(linspace(10,100,ngrid)', linspace(0,H*maxshare,ngrid)');

price					=	zeros(ngrid,ngrid);

for i=1:ngrid
	if mod(i,10)==0
		fprintf('Loop %d of %d\n',i,ngrid);
	end
	for j=1:ngrid
		price(i,j)	=	norminv(1-(H - invgrid(i,j))/N,0,bmgrid(i,j));
	end
end

f	=	figure;
%  set(f,'Renderer','painters');
hs	=	surf(bmgrid,invgrid,price,'FaceAlpha',0.1);
alpha(0.1);
%  h1	=	contour(bmgrid,invgrid,price,'Parent',gca);
%  alpha(0.1);
xlabel('Bid Variance');
ylabel('Inventory');
zlabel('Price');
title('Variance/inventory analysis');
print(f,'-djpeg','outgraph/biddist.jpg');
%  close(f);

