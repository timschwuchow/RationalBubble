%  Copyright 2012 Timothy John Schwuchow
%  psimpolgraph002.m		-	Graph policy functions

load('outmat/specparam.mat');

ieval	=	1;
eeval	=	ceil(epts/2);
deval	=	ceil(dpts/2);
evec	=	(1:epts);
dvec	=	(1:dpts);
uvec	=	(1:upts);
for i=1:ngrid
	ii=i;
	fprintf('Graphing trial %d \n \n', i);
	load(['outmat/pol' int2str(i) '.mat']);
	i=ii;
	esurf	=	zeros(epts,1);
	usurf	=	zeros(upts,1);
	isurf	=	zeros(upts,epts);
	for j=1:epts
		esurf(j)	=	emap(i,j);
		for k=1:upts
			usurf(k)	=	umap(i,k);
			isurf(k,j)	=	imap(pol(ieval,smap(j,k,deval)));
		end
		if (mod(j,5)==0)
			fprintf('UE graph %4.2f %% percent complete \n', 100*j/epts);
		end
	end
	f	=	figure;
	set(gcf, 'PaperPositionMode','manual');
	set(gcf, 'PaperUnits', 'inches');
	set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
	h	=	surfc(esurf,usurf,isurf);
	xlabel('\eta');
	ylabel('u');
	zlabel('h');
	% print(f,'-djpeg',['out1/pol' int2str(i) 'ue.jpg']);
	close(f);
	dsurf	=	zeros(dpts,1);
	usurf	=	zeros(upts,1);
	isurf	=	zeros(upts,dpts);
	for j=1:dpts
		dsurf(j)	=	dmap(i,j);
		for k=1:upts
			usurf(k)	=	umap(i,k);
			isurf(k,j)	=	imap(pol(ieval,smap(eeval,k,j)));
		end
		if mod(j,5)==0
			fprintf('UD graph %4.2f %% complete \n', 100*j/dpts);
		end
	end
	f	=	figure;
	set(gcf, 'PaperPositionMode','manual');
	set(gcf, 'PaperUnits', 'inches') ;
	set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
	h	=	surfc(dsurf,usurf,isurf);
	xlabel('\delta');
	ylabel('u');
	zlabel('h');
	print(f,'-djpeg',['out1/pol' int2str(i) 'ud.jpg']);
	close(f);
end
