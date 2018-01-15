%  Copyright 2012 Timothy John Schwuchow
%  psimspecstats002.m		-	Statistics for consumer-only and speculator models:

clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');

load('outmat/specparam.mat');
pspecall			=	zeros(nsim,T,numel(regrid));
pconsall			=	pspecall;
etaall			=	pspecall;
deltacons		=	pspecall;
innoall			=	pspecall;
deltaspec		=	pspecall;
invall			=	pspecall;
cfall				=	pspecall;
piall				=	pspecall;
ngrid				=	numel(regrid);
corrlags			=	16;

disp('Load')
for k=1:ngrid
	fprintf('Loading simulation %d \n',k);
	load(['outmat/pricesim' int2str(k) '.mat']);
	pspecall(:,:,k)		=	pspec;
	pconsall(:,:,k)		=	pcons;
	peffall(:,:,k)			=	peff;
%  	etaall(:,:,i)			=	reshape(emap(smapinv(statespec,1)),size(statespec));
%  	innoall(:,:,i)			=	reshape(umap(smapinv(statespec,2)),size(statespec));
%  	deltacons(:,:,i)		=	reshape(dmap(smapinv(statecons,3)),size(statespec));
%  	deltaspec(:,:,i)		=	reshape(dmap(smapinv(statespec,3)),size(statespec));
%  	max(inv(1,:))
	invall(:,:,k)			=	inv;
	cfall(:,:,k)			=	cf;
	piall(:,:,k)			=	pi;
end

%  Create Basic statistics and graphs
disp('Generate Tables')
pvec				=	'Simulation & $\sigma_\eta$ & $\sigma_\gamma$ & $\rho_\eta$ & $\rho_\gamma$';
outtab11			=	fopen(['out1/invpatterns.tex'], 'w');
outtab12			=	fopen(['out1/profits.tex'], 'w');
outtab13			=	fopen(['out1/momentum.tex'], 'w');
outtab14			=	fopen(['out1/priceeff.tex'], 'w');

%  Headers for tables
fprintf(outtab11,['\\begin{tabular}{c | c c c c | c c | c c}\n %s ' repmat(' & %s ',1,4) ' \\\\ \\hline \n'],pvec,'$\bar{h}_t$','$Var(h_t)$','$\bar{\Delta h_t}$','$Var(\Delta h_t)$');

fprintf(outtab12,['\\begin{tabular}{c | c c c c | c c | c c}\n %s ' repmat(' & %s ',1,4) ' \\\\ \\hline \n'],pvec,'$\bar{R}_t$','$Var(R_t)$','$\bar{CF_t}$','$Var(CF_t)$');
corrvec			=	[];
for i=1:corrlags
	corrvec		=	[corrvec ' & $\hat{\rho}_{' int2str(i) '}$ '];
end
fprintf(outtab13,['\\begin{tabular}{c c | c c c c | ' repmat(' c ',1,corrlags) ' }\n %s & %s %s \\\\ \\hline \n'],'Model',pvec,corrvec);


fprintf(outtab14,['\\begin{tabular}{c | c c c c | ' repmat(' c ',1,6) ' }\n %s ' repmat(' & %s ',1,6) ' \\\\ \\hline \n'],pvec,'$\frac{\bar{p}_c - \bar{p}_e}{\bar{p}_e}$','$\frac{\bar{p}_s - \bar{p}_e}{\bar{p}_e}$','$E\left[abs\left(\frac{p_c - p_e}{p_e}\right) \right]$','$E\left[abs\left(\frac{p_s - p_e}{p_e}\right) \right]$','$E\left[\frac{Var(\Delta p_c) - Var(\Delta p_e)}{Var(\Delta p_e)}\right]$','$E\left[\frac{Var(\Delta p_s) - Var(\Delta p_e)}{Var(\Delta p_e)}\right]$');

for i=1:ngrid
	fprintf('Generate stats for simulation %d \n',i);
	paramlabel		=	sprintf('$\\sigma_\\eta$ = %5.3f | $\\sigma_\\gamma$ = %5.3f | $\\rho_\\eta$ = %5.3f | $\\rho_\\gamma$ = %5.3f', pgrid(i,:));
	% Inventory
	longinv			=	reshape(invall(:,:,i),[numel(inv) 1]);
	meaninv			=	mean(longinv);
	varinv			=	var(longinv);
	% Profits and Cash Flow
	rofr				=	log(1+piall(:,2:T,i)-piall(:,1:T-1,i));
	meanrofr			=	mean(reshape(rofr,[numel(rofr) 1]));
	varrofr			=	var(reshape(rofr,[numel(rofr) 1]));
	longcf			=	reshape(cfall(:,:,i),[numel(cfall(:,:,i)) 1]);
	meancf			=	mean(longcf);
	varcf				=	var(longcf);
	[invdensityy invdensityx]	=	ksdensity(longinv);
	% Inventory Changes
	deltainv			=	invall(:,2:T,i) - invall(:,1:T-1,i);
	longdeltainv	=	reshape(deltainv,[numel(deltainv) 1]);
	meandeltainv	=	mean(longdeltainv);
	vardeltainv		=	var(longdeltainv);
	[deltainvdensityy deltainvdensityx]	=	ksdensity(longdeltainv);

	% Price differences
	pspecdif			=	pspecall(:,2:T,i) - pspecall(:,1:T-1,i);
	pconsdif			=	pconsall(:,2:T,i) - pconsall(:,1:T-1,i);
	peffdif			=	peffall(:,2:T,i) - peffall(:,1:T-1,i);
	longpspecdif	=	reshape(pspecdif,[nsim*(T-1) 1]);
	longpconsdif	=	reshape(pconsdif,[nsim*(T-1) 1]);
	longpeffdif		=	reshape(peffdif,[nsim*(T-1) 1]);

	% Price Difference momentum
	corrpspec		=	zeros(1,corrlags);
	corrpcons		=	zeros(1,corrlags);
	corrpeff			=	zeros(1,corrlags);

	for j=1:corrlags
		cf					=	corrcoef(pspecdif(:,1+j:T-1),pspecdif(:,1:T-1-j));
		corrpspec(j)	=	cf(2);
		cf					=	corrcoef(pconsdif(:,1+j:T-1),pconsdif(:,1:T-1-j));
		corrpcons(j)	=	cf(2);
		cf					=	corrcoef(peffdif(:,1+j:T-1),peffdif(:,1:T-1-j));
		corrpeff(j)		=	cf(2);
	end
	% Prices relative to efficient
	longpdifce		=	reshape((pconsall(:,:,i) - peffall(:,:,i))./abs(peffall(:,:,i)),[nsim*T 1]);
	longpdifse		=	reshape((pspecall(:,:,i) - peffall(:,:,i))./abs(peffall(:,:,i)),[nsim*T 1]);
	meanpdifce		=	mean(longpdifce);
	meanpdifse		=	mean(longpdifse);
	ameanpdifce		=	mean(abs(longpdifce));
	ameanpdifse		=	mean(abs(longpdifse));

	% Relative volatility (variance of price changes)
	varpdifse		=	(var(reshape(pspecall(:,:,i),[nsim*T 1])) - var(reshape(peffall(:,:,i),[nsim*T 1])))/var(reshape(peffall(:,:,i),[nsim*T 1]));
	varpdifce		=	(var(reshape(pconsall(:,:,i),[nsim*T 1])) - var(reshape(peffall(:,:,i),[nsim*T 1])))/var(reshape(peffall(:,:,i),[nsim*T 1]));



	if i==ngrid
		endlines		=	'\\hline';
	else
		endlines		=	'';
	end
	fprintf(outtab11,['%d ' repmat(' & %5.3f ', 1,8) '\\\\' endlines endlines ' \n'],i,pgrid(i,:),meaninv,varinv,meandeltainv,vardeltainv);
	fprintf(outtab12,['%d ' repmat(' & %5.3f ', 1,8) '\\\\' endlines endlines ' \n'],i,pgrid(i,:),meanrofr,varrofr,meancf,varcf);

%  	fprintf(outtab13,['%s ' repmat(' & %5.3f ', 1,corrlags+4) '\\\\ \n'],'Efficient',pgrid(i,:),corrpeff);
	fprintf(outtab13,['%s & %d ' repmat(' & %5.3f ', 1,corrlags+4) '\\\\ \n'],'No Spec.',i,pgrid(i,:),corrpcons);
	fprintf(outtab13,['%s & %d ' repmat(' & %5.3f ', 1,corrlags+4) '\\\\ \\hline ' endlines ' \n'],'Speculator',i,pgrid(i,:),corrpspec);

	fprintf(outtab14,['%d & %5.3f ' repmat(' & %5.3f ',1,9) '\\\\' endlines endlines ' \n'],i,pgrid(i,:),meanpdifce,meanpdifse,ameanpdifce,ameanpdifse,varpdifce,varpdifse);
	f1	=	figure;
	set(gcf, 'PaperPositionMode','manual');
	set(gcf, 'PaperUnits', 'inches') ;
	set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
	h1	=	plot(invdensityx,invdensityy,'-r');
	xlabel('Inventory');
	ylabel('Density');
	title(['Inventory Holdings:' paramlabel]);
	print(f1,'-djpeg',['out1/invdensity' int2str(i) '.jpg']);
	close(f1);

	f1	=	figure;
	set(gcf, 'PaperPositionMode','manual');
	set(gcf, 'PaperUnits', 'inches') ;
	set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
	h1	=	plot(deltainvdensityx,deltainvdensityy,'-r');
	xlabel('Inventory Changes');
	ylabel('Density');
	title(['Inventory Changes:' paramlabel]);
	print(f1,'-djpeg',['out1/deltainvdensity' int2str(i) '.jpg']);
	close(f1);
end
fprintf(outtab11,'\\end{tabular} \n');
fclose(outtab11);
fprintf(outtab12,'\\end{tabular} \n');
fclose(outtab12);
fprintf(outtab13,'\\end{tabular} \n');
fclose(outtab13);
fprintf(outtab14,'\\end{tabular} \n');
fclose(outtab14);