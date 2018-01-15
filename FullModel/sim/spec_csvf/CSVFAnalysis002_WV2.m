%  Copyright 2012-2013 Timothy John Schwuchow
%  CSVFAnalysis002_WV2.m		-	Analyze CSVF variation with parameters (start from base case)
%  Working version
clear all
format short g
path(path,'../../functions');
%  Control parameters
seed		=	8;
T			=	1000;
nsim		=	100;
N			=	1;
epts		=	65;
dpts		=	65;
upts		=	65;
ipts		=	40;
esd			=	4.00;
usd			=	3.00;
dsd			=	3.00;
maxshare	=	0.90;
vftolcrit	=	5e-5;
vfmaxiter	=	5000;
howimp		=	40;
howlag		=	10;
nonaive		=	0;

%  Simulation parameters
H		=	[0.500];
x		=	[00.00];
xvar	=	[0.000];
ve		=	[1.000];
vg		=	[10.00];
re		=	[0.990];
rg		=	[0.250];
b		=	0.99;
vdiag	=	0;
paramcurrent	=	[ve,vg,re,rg,x,xvar,H,b];
fprintf('Initializing CSVF Interactive Analyzer\nStates: %d\n',ipts*dpts*upts*epts);
passloop = 0;
ct	=	0;
if exist('outmat/csvfbase.mat')==2
	if exist('outcsv/parameters.csv')==2
		pold	=	dlmread('outcsv/parameters.csv','\t',1,0);
		if sum([paramcurrent epts upts dpts ipts]~=pold(1,2:numel(paramcurrent)+5))>0
			delete('outmat/csvfbase.mat');
			delete('outcsv/parameters.csv');
			CSVFBase				=	cell(ipts,1);
			VBase					=	cell(ipts,1);
			PolBase				=	cell(ipts,1);
			for i=1:ipts
				VBase{i}			=	zeros(1,epts*upts*dpts);
				PolBase{i}		=	zeros(1,epts*upts*dpts);
				CSVFBase{i}		=	zeros(ipts,epts*upts*dpts);
			end
		else
			fprintf('Loading base matrix\n');
			load outmat/csvfbase.mat
			CSVFBase	=	CSVFCurrent;
			VBase		=	VCurrent;
			PolBase	=	PolCurrent;
			passloop = 	1;
			ct			=	size(pold,1);
		end
	else
		CSVFBase				=	cell(ipts,1);
		VBase					=	cell(ipts,1);
		PolBase				=	cell(ipts,1);
		for i=1:ipts
			VBase{i}			=	zeros(1,epts*upts*dpts);
			PolBase{i}		=	zeros(1,epts*upts*dpts);
			CSVFBase{i}		=	zeros(ipts,epts*upts*dpts);
		end
	end
else
	if exist('outcsv/parameters.csv')==2
		delete('outcsv/parameters.csv');
	end
	CSVFBase				=	cell(ipts,1);
	VBase					=	cell(ipts,1);
	PolBase				=	cell(ipts,1);
	for i=1:ipts
		VBase{i}			=	zeros(1,epts*upts*dpts);
		PolBase{i}		=	zeros(1,epts*upts*dpts);
		CSVFBase{i}		=	zeros(ipts,epts*upts*dpts);
	end
end
fprintf('Finished Initializing matrices\n')




loopind			=	'Y';
vars			=	{'ve','vg','re','rg','x','xvar','H','b'};
msg				=	{'Set value for '};
msg				=	strcat(msg, vars,{': [Default '});
parambase		=	paramcurrent;
fprintf('Begin CSVF Loop\n');
while (loopind == 'Y')
	if passloop ~= 1
		tloop						=	tic;
		% Initialize
		ct							=	ct+1;
		paramcurrent			=	zeros(1,8);
		CSVFCurrent				=	CSVFBase;
		VCurrent					=	VBase;
		PolCurrent				=	PolBase;
		fprintf('Parameter history\n');
		unix('cat outcsv/parameters.csv; echo ""');
		for i=1:8
			if ct ~=1
				ntemp					=	input([msg{i} num2str(parambase(i),'%4.2f') '] '], 's');
				if not(isempty(ntemp))
					paramcurrent(i)	=	str2num(ntemp);
				else
					paramcurrent(i)	=	parambase(i);
				end
			else
				paramcurrent(i)	=	parambase(i);
			end
		end
		ve		=	paramcurrent(1);
		vg		=	paramcurrent(2);
		re		=	paramcurrent(3);
		rg		=	paramcurrent(4);
		x		=	paramcurrent(5);
		xvar	=	paramcurrent(6);
		H		=	paramcurrent(7);
		b		=	paramcurrent(8);

		[cd1,ce1,cd2,ce2,cb2,bmvar,imap,emap,dmap,umap,smap,smapinv]	=	CSVFParamSet002(ve,vg,re,rg,x,xvar,b,vdiag,H,N,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,ct);
		[CSVFCurrent VCurrent PolCurrent]	=	 CSVFUpdate002(CSVFBase,VBase,PolBase,ve,re,cd1,ce1,cd2,ce2,cb2,bmvar,b,H,N,x,vftolcrit,vfmaxiter,imap,emap,dmap,umap,smap,smapinv,howimp,howlag,nonaive);

		csvfsurface1	=	zeros(ipts,upts);
		csvfsurface2	=	zeros(ipts,dpts);
		[~,dzero]		=	min(abs(0 - dmap));
		[~,uzero]		=	min(abs(0 - umap));
		for yy=1:ipts
			for xx=1:upts
				csvfsurface1(yy,xx)	=	CSVFCurrent{1}(yy,smap(ceil(epts/2),xx,dzero));
			end
			for xx=1:dpts
				csvfsurface2(yy,xx)	=	CSVFCurrent{1}(yy,smap(ceil(epts/2),uzero,xx));
			end
		end
		maxval1			=	zeros(upts,1);
		invopt1			=	zeros(upts,1);
		for xx=1:upts
			invopt1(xx)	=	imap(PolCurrent{1}(smap(ceil(epts/2),xx,dzero)));
			maxval1(xx)	=	VCurrent{1}(smap(ceil(epts/2),xx,dzero));
		end
		maxval2			=	zeros(dpts,1);
		invopt2			=	zeros(dpts,1);
		for xx=1:dpts
			invopt2(xx)	=	imap(PolCurrent{1}(smap(ceil(epts/2),uzero,xx)));
			maxval2(xx)	=	VCurrent{1}(smap(ceil(epts/2),uzero,xx));
		end


		f1		=	figure;
		h1		=	surfc(umap,imap,csvfsurface1);
		ax1	=	gca;
		hold on
		h2		=	plot3(umap,invopt1,maxval1,'-k','Parent',ax1,'LineWidth',2.0);
		hold on
		ZRange	=	get(gca,'ZLim');
		h3		=	plot3(umap,invopt1,ZRange(1)*ones(size(umap)),'-k','Parent',ax1,'LineWidth',2.0);
		title(sprintf('Trial %d :: \\mu CSVF :: \\sigma_\\mu^2 = %5.2f, \\sigma_\\xi^2 = %5.2f, \\rho_\\eta = %5.2f , \\rho_\\gamma = %5.2f, x = %5.2f, \\sigma_x^2 = %5.2f, H = %5.2f, b = %5.2f',ct,ve,vg,re,rg,x,xvar,H,b));
		xlabel('\mu');
		ylabel('Inv.');
		zlabel('Value');
		print(f1,'-djpeg',['outgraph/CSVFMu' num2str(ct,'%d') '.jpg']);
		close(f1);
		f1	=	figure;
		set(f1,'Renderer','painters');
		h1	=	surfc(dmap,imap,csvfsurface2);
		ax1	=	gca;
		hold on
		h2		=	plot3(dmap,invopt2,maxval2,'-k','Parent',ax1,'LineWidth',2.0);
		hold on
		ZRange	=	get(gca,'ZLim');
		h3		=	plot3(dmap,invopt2,ZRange(1)*ones(size(dmap)),'-k','Parent',ax1,'LineWidth',2.0);
		title(sprintf('Trial %d :: \\delta CSVF :: \\sigma_\\mu^2 = %5.2f, \\sigma_\\xi^2 = %5.2f, \\rho_\\eta = %5.2f , \\rho_\\gamma = %5.2f, x = %5.2f, \\sigma_x^2 = %5.2f, H = %5.2f, b = %5.2f',ct,ve,vg,re,rg,x,xvar,H,b));
		xlabel('\delta');
		ylabel('Inv.');
		zlabel('Value');
		print(f1,'-djpeg',['outgraph/CSVFDelta' num2str(ct,'%d') '.jpg']);
		close(f1);
		fprintf('Loop Completed in %4.2f minutes\n',toc(tloop)/60);
		loopind	=	upper(input('Continue Interactive Analyzer? Y/N [Y]: ', 's'));
		if isempty(loopind) | loopind ~= 'N'
			loopind	=	'Y';
		end
		if ct==1 & exist('outmat/csvfbase.mat')~=2
			save('outmat/csvfbase.mat','CSVFCurrent','VCurrent','PolCurrent','-V7.3');
		end
		passloop=0;
	else
		fprintf('Passing first loop\n');
		passloop=0;
	end
end
%

