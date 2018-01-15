%  Copyright 2012-2013 Timothy John Schwuchow
%  CSVFAnalysis002_WV.m		-	Analyze CSVF variation with parameters
%  Working version
clear all
format short g
path(path,'../../functions');
%  Control parameters
seed		 	=	8;
T				=	1000;
nsim			=	100;
N				=	1;
epts			=	51;
dpts			=	51;
upts			=	41;
ipts			=	40;
esd			=	4.00;
usd			=	3.00;
dsd			=	3.00;
maxshare		=	0.90;
vftolcrit	=	1e-4;
vfmaxiter	=	5000;
howimp		=	10;
howlag		=	10;
nonaive		=	0;

%  Simulation parameters
H		=	0.50;
x		=	[10.00];
xvar	=	[0.000];
ve		=	[1.000];
vg		=	[1.050];
re		=	[0.990];
rg		=	[0.250];
b		=	0.99;
vdiag	=	0;

fprintf('Initializing CSVF Interactive Analyzer\n');
if exist('outmat/csvfinit.mat')==2
	load outmat/csvfinit.mat
else
	CSVFCurrent				=	cell(ipts,1);
	VCurrent					=	cell(ipts,1);
	PolCurrent				=	cell(ipts,1);
	for i=1:ipts
		VCurrent{i}			=	zeros(1,epts*upts*dpts);
		PolCurrent{i}		=	zeros(1,epts*upts*dpts);
		CSVFCurrent{i}		=	zeros(ipts,epts*upts*dpts);
	end
end
fprintf('Finished Initializing matrices\n')



if exist('outcsv/parameters.csv')==2
	delete('outcsv/parameters.csv');
end

loopind			=	'Y';
vars				=	{'ve','vg','re','rg','x','xvar','H','b'};
msg				=	{'Set value for '};
msg				=	strcat(msg, vars,{': [Default '});
paramcurrent	=	[ve,vg,re,rg,x,xvar,H,b];
ct					=	0;
fprintf('Begin CSVF Loop\n');
while (loopind == 'Y')
	% Initialize
	ct							=	ct+1;
	paramlast				=	paramcurrent;
	paramcurrent			=	zeros(8,1);
	CSVFLast					=	CSVFCurrent;
	VLast						=	VCurrent;
	PolLast					=	PolCurrent;
	for i=1:8
		if ct ~=1
			ntemp					=	input([msg{i} num2str(paramlast(i),'%4.2f') '] '], 's');
			if not(isempty(ntemp))
				paramcurrent(i)	=	str2num(ntemp);
			else
				paramcurrent(i)	=	paramlast(i);
			end
		else
			paramcurrent(i)	=	paramlast(i);
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

	[cd1,ce1,cd2,ce2,cb2,bmvar,imap,emap,dmap,umap,smap,smapinv]	=	CSVFParamSet002(ve,vg,re,rg,x,xvar,b,vdiag,H,N,ipts,epts,upts,dpts,maxshare,esd,usd,dsd);
	[CSVFCurrent VCurrent PolCurrent]	=	 CSVFUpdate002(CSVFLast,VLast,PolLast,ve,re,cd1,ce1,cd2,ce2,cb2,bmvar,b,H,N,x,vftolcrit,vfmaxiter,imap,emap,dmap,umap,smap,smapinv,howimp,howlag,nonaive);

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
%  	set(f1,'Renderer','painters');
%  	ax1	=	axes('XLim',[min(umap) max(umap)],'YLim',[min(imap) max(imap)],'ZLim',[m max(max(csvfsurface1))], 'CLim',[min(min(csvfsurface1)) max(max(csvfsurface1))]);
	h1		=	surfc(umap,imap,csvfsurface1);
	ax1	=	gca;
	hold on
	h2		=	plot3(umap,invopt1,maxval1,'-k','Parent',ax1,'LineWidth',2.0);
	hold on
	ZRange	=	get(gca,'ZLim');
	h3		=	plot3(umap,invopt1,ZRange(1)*ones(size(umap)),'-k','Parent',ax1,'LineWidth',2.0);
	titlestring = sprintf('Trial %d :: \\mu CSVF for \\sigma_\\mu^2 = %4.2f, \\sigma_\\xi^2 = %4.2f, \\rho_\\eta = %4.2f , \\rho_\\gamma = %4.2f, x = %4.2f, \\sigma_x^2 = %4.2f, N = %4.2f, b = %4.2f',ct,ve,vg,re,rg,x,xvar,N,b);
	title(titlestring);
	xlabel('\mu');
	ylabel('Inv.');
	zlabel('Value');
	print(f1,'-djpeg',['outgraph/CSVFMu' num2str(ct,'%d') '.jpg']);
	close(f1);
	f1	=	figure;
%  	set(f1,'Renderer','painters');
	h1	=	surfc(dmap,imap,csvfsurface2);
	ax1	=	gca;
	hold on
	h2		=	plot3(dmap,invopt2,maxval2,'-k','Parent',ax1,'LineWidth',2.0);
	hold on
	ZRange	=	get(gca,'ZLim');
	h3		=	plot3(dmap,invopt2,ZRange(1)*ones(size(dmap)),'-k','Parent',ax1,'LineWidth',2.0);
	titlestring = sprintf('Trial %d :: \\delta CSVF for \\sigma_\\mu^2 = %5.2f, \\sigma_\\xi^2 = %5.2f, \\rho_\\eta = %5.2f , \\rho_\\gamma = %5.2f, x = %5.2f, \\sigma_x^2 = %5.2f, N = %5.2f, b = %5.2f',ct,ve,vg,re,rg,x,xvar,N,b));
	title(titlestring);
	xlabel('\delta');
	ylabel('Inv.');
	zlabel('Value');
	print(f1,'-djpeg',['outgraph/CSVFDelta' num2str(ct,'%d') '.jpg']);
	close(f1);
	loopind	=	upper(input('Continue Interactive Analyzer? Y/N [Y]: ', 's'));
	if isempty(loopind)
		loopind	=	'Y';
	end
	if ct==1 & exist('outmat/csvfinit.mat')~=2
		save('outmat/csvfinit.mat','CSVFCurrent','VCurrent','PolCurrent','-V7.3');
	end
end
%

