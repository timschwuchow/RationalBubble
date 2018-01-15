%  Copyright 2012 Timothy John Schwuchow
clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');

load('../../sim/t4/outmat/pricesim1.mat');
trial			=	1;

blarg			= 	zeros(epts,3);
elong			=	reshape(emap(smapinv(statecons,1)),numel(statecons),1);
for i=1:epts
	blarg(i,:)	=	[emap(i) sum(elong==emap(i)) sum(elong==emap(i))/numel(elong)];
end
blarg;
bmcons	=	norminv(1 - H/ N,0,bmvargrid);
d_T			=	T-1;
dd_T			=	T-2;

pconshat		=	(1-b*regrid)*(pcons(trial,:)' - (x + bmcons)/(1-b));
d_pconshat	=	pconshat(2:T) - regrid*pconshat(1:T-1);
delta			=	dmap(smapinv(statecons(trial,:),3))';
errs			=	umap(smapinv(statecons(trial,:),2))';
eta			=	emap(smapinv(statecons(trial,:),1))';
d_errcfs		=	[1 + b*regrid*ce1grid, b*regrid*cd1grid*ce2grid - b*regrid*ce1grid*regrid, b*regrid*cd1grid*cd2grid - b*regrid*cd1grid*regrid];
d_errlags	=	[errs(2:T), errs(1:T-1), delta(1:T-1)];
compare		=	[d_pconshat sum(repmat(d_errcfs,T-1,1).*d_errlags,2)];
deltaerr		=	delta(2:T) - cd2grid*delta(1:T-1) - ce2grid*errs(1:T-1);
etaerr		=	eta(2:T) - regrid*eta(1:T-1) - errs(2:T);
[compare(1:10,:) compare(1:10,1)-compare(1:10,2)];
for i=1:epts
	blarg(i,:)	=	[emap(i) sum(eta==emap(i)) sum(eta==emap(i))/numel(eta)];
end
blarg;


dd_pconshat	=	d_pconshat(2:T-1) - cd2grid*d_pconshat(1:d_T-1);
errcfs		=	[ (1 + b*regrid*ce1grid) , -cd2grid - b*regrid*cd2grid*ce1grid + b*regrid*cd1grid*ce2grid - b*regrid^2*ce1grid, b*regrid^2*(ce1grid*cd2grid - cd1grid*ce2grid)];
ulags			=	[errs(3:T) errs(2:T-1) errs(1:T-2)];
rhs			=	sum(repmat(errcfs,T-2,1).*ulags,2);
[dd_pconshat(1:10),rhs(1:10),dd_pconshat(1:10)-rhs(1:10)];
dif			=	dd_pconshat - rhs;

TT				=	T-2;
samp			=	mean(dd_pconshat);
for i=0:5
	cf			=	cov(dd_pconshat(1+i:TT),dd_pconshat(1:TT-i));
	samp		=	[samp,cf(2)];
end
samp
%  samp=[mean(dd_pconshat) var(dd_pconshat) cov(dd_pconshat(2:TT),dd_pconshat(1:TT-1)),cov(dd_pconshat(3:TT),dd_pconshat(1:TT-2)),cov(dd_pconshat(4:TT),dd_pconshat(1:TT-3)),cov(dd_pconshat(5:TT),dd_pconshat(1:TT-4))]

%  dif(1:10)
fprintf('Single Trial:: Mean dif = %7.5g | RMSE = %7.5g | Var = %7.5g \n',mean(dif),mean(dif.^2)^0.5, var(dif))

clear all

load('../../sim/t4/outmat/pricesim1.mat');

bmcons	=	norminv(1 - H/ N,0,bmvargrid);


pconshat		=	(1-b*regrid)*(pcons - (x + bmcons)/(1-b));
dd_pconshat	=	reshape(pconshat(:,3:T)' - regrid*pconshat(:,2:T-1)' - cd2grid*(pconshat(:,2:T-1)' - regrid*pconshat(:,1:T-2)'),[numel(pconshat(:,3:T)) 1]);



errs			=	umap(smapinv(statecons',2))';
errcfs		=	[ (1 + b*regrid*ce1grid) , -cd2grid - b*regrid*cd2grid*ce1grid + b*regrid*cd1grid*ce2grid - b*regrid^2*ce1grid, b*regrid^2*(ce1grid*cd2grid - cd1grid*ce2grid)];
ulags			=	[umap(smapinv(statecons(:,3:T)',2))',umap(smapinv(statecons(:,2:T-1)',2))',umap(smapinv(statecons(:,1:T-2)',2))'];
[dd_pconshat(1:10) ulags(1:10,:)];
rhs			=	sum(repmat(errcfs,size(ulags,1),1).*ulags,2);

dif			=	dd_pconshat - rhs;
%  dif(1:10)
fprintf('All Trials:: Mean dif = %7.5g | RMSE = %7.5g | Var = %7.5g \n',mean(dif),mean(dif.^2)^0.5, var(dif))

errcfs		=	[ (1 + b*regrid*ce1grid) , -cd2grid - b*regrid*cd2grid*ce1grid + b*regrid*cd1grid*ce2grid - b*regrid^2*ce1grid, b*regrid^2*(ce1grid*cd2grid - cd1grid*ce2grid)];

moms			=	5;
errcfs		=	[errcfs, zeros(1,moms)];
tmoms			=	[0];
%  tmoms			=	[0,sum(errcfs(1:3).^2).vegrid];
for i=0:moms
	tmoms		=	[tmoms sum(vegrid*errcfs(1:3).*errcfs(1+i:3+i))];
end

true = tmoms
