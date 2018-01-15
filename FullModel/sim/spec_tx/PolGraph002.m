%  Copyright 2012-2013 Timothy John Schwuchow
%  PolGraph002.m		-	Graph Policy Function
%  Production version

clear all

format short g
path(path,'../../functions');
load('outmat/specparam.mat');
load(['outmat/pol' int2str(repng) '.mat']);

uvec	=	(1:upts)';
dvec	=	(1:dpts)';

polsurface	=	zeros(dpts,upts);

for j=1:upts
	for i=1:dpts
		polsurface(i,j)	=	imap(pol(1,smap(ceil(epts/2),j,i)));
	end
end

f1	=	figure;
set(f1,'Renderer','painters')
h1	=	surfc(umap(repng,:),dmap(repng,:),polsurface);
title(sprintf('Policy Function for \\sigma_\\mu^2 = %5.3f, \\sigma_\\xi^2 = %5.3f, \\rho_\\eta = %5.3f , \\rho_\\gamma = %5.3f, x = %5.3f, \\sigma_x^2 = %5.3f',pgrid(repng,:)));
xlabel('\mu');
ylabel('\delta');
zlabel('h');
print(f1,'-djpeg','outgraph/policyrepng.jpg')
quit
