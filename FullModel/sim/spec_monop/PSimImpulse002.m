%  Copyright 2012-2013 Timothy John Schwuchow
%  PSimImpulse002.m		-	Generate impulse response and graph
%  Production version

clear all

format short g
path(path,'../../functions');
load('outmat/specparam.mat');
load(['outmat/pol' int2str(repng) '.mat']);

randn('state',seed);
pgrid			=	pgrid(repng,:);
emap			=	emap(repng,:);
dmap			=	dmap(repng,:);
umap			=	umap(repng,:);
regrid		=	regrid(repng);
x				=	xgrid(repng);
bmvargrid	=	bmvargrid(repng);
varvgrid		=	varvgrid(repng);
cd1grid		=	cd1grid(repng);
ce1grid		=	ce1grid(repng);
vegrid		=	vegrid(repng);
cd2grid		=	cd2grid(repng);
ce2grid		=	ce2grid(repng);
cb2grid		=	cb2grid(repng);
tg				=	30;
split			=	ceil(tg/2);
%  Speculative Prices

pspec				=	zeros(2,tg);
statespec		=	zeros(2,tg);
[~,uinit]		=	min(abs(umap - 0));
[~,dinit]		=	min(abs(dmap - 0));
[~,einit]		=	min(abs(emap - 0));
statespec(:,1)	=	smap(einit,uinit,dinit);
invstate			=	zeros(2,tg);
inv				=	zeros(2,tg);
iinit				=	[1;5];

inno				=	uinit*ones(2,tg);
inno(1,split)	=	[ceil(0.7*upts)];
inno(2,split)	=	[ceil(0.3*upts)];
for t=1:tg
	if t<split
		invstate(:,t)			=	iinit;
	else
		invstate(:,t)			=	pol(sub2ind(size(pol),invstate(:,t-1),statespec(:,t)));

	end
	inv(:,t)						=	imap(invstate(:,t))';
	pspec(:,t)					=	pricespec002(inv(:,t),[0;0],emap(smapinv(statespec(:,t),1))',dmap(smapinv(statespec(:,t),3))',umap(smapinv(statespec(:,t),2))',x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
	if t < tg
		for i=1:2
			naiveerr				=	(norminv(1 - (H)/ N,0,bmvargrid) - norminv(1 - (H-inv(i,t))/ N,0,bmvargrid));
			[~,etemp]			=	min(abs(emap - (regrid*emap(smapinv(statespec(i,t),1))' + umap(inno(i,t+1)))),[],2);
			[~,dtemp]			=	min(abs(dmap - (cb2grid*naiveerr + cd2grid*dmap(smapinv(statespec(i,t),3)) + ce2grid*umap(smapinv(statespec(i,t),2)))),[],2);
			statespec(i,t+1)	=	smap(etemp,inno(i,t+1),dtemp);
		end
	end
end


% Competitive Speculator Prices
PriceSpecComp			=	zeros(2,tg);
StateSpecComp			=	zeros(2,tg);
StateSpecComp(:,1)	=	statespec(:,1);
SpecInvComp				=	zeros(2,tg);
SpecInvComp(:,1)		=	inv(:,1);
thresh					=	1-1e-10;
for t=1:tg
	if t > 1
		for i=1:2
			utemp						=	smapinv(statespec(i,t),2);
			etemp						=	smapinv(statespec(i,t),1);
			[~,dtemp]				=	min(abs(dmap - (cd2grid*dmap(smapinv(StateSpecComp(i,t-1),3)) + ce2grid*umap(smapinv(StateSpecComp(i,t-1),2)) + cb2grid*(norminv(1-H/N,0,bmvargrid) - norminv(1-(H-SpecInvComp(i,t-1))/N,0,bmvargrid)))),[],2);
			StateSpecComp(i,t)	=	smap(etemp,utemp,dtemp);
		end
		if t < split
			SpecInvComp(:,t)		=	SpecInvComp(:,t-1);
		else
			SpecInvComp(:,t)		=	max(N*(normcdf( -(x+b*norminv(1-H/N,0,bmvargrid))/(1-b) - emap(smapinv(StateSpecComp(:,t),1))'/(1-b*regrid) - cd1grid*dmap(smapinv(StateSpecComp(:,t),3))' - ce1grid*umap(smapinv(StateSpecComp(:,t),2))',0,bmvargrid)-1)+H*maxshare,0);
		end
	end
	PriceSpecComp(:,t)		=	pricespec002(SpecInvComp(:,t),zeros(2,1),emap(smapinv(StateSpecComp(:,t),1))',dmap(smapinv(StateSpecComp(:,t),3))',umap(smapinv(StateSpecComp(:,t),2))',x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
end

%  Consumer Prices

statecons				=	zeros(2,tg);
statecons(:,1)			=	statespec(:,1);
pcons						=	zeros(2,tg);
for t=1:tg
	pcons(:,t)			=	pricespec002([0;0],[0;0],emap(smapinv(statecons(:,t),1))',dmap(smapinv(statecons(:,t),3))',umap(smapinv(statecons(:,t),2))',x,b,regrid,H,N,bmvargrid,cd1grid,ce1grid);
	if t < tg
		for i=1:2
			utemp					=	smapinv(statespec(i,t+1),2);
			etemp					=	smapinv(statespec(i,t+1),1);
			[~,dtemp]			=	min(abs(dmap - (cd2grid*dmap(smapinv(statecons(i,t),3)) + ce2grid*umap(smapinv(statecons(i,t),2)))),[],2);
			statecons(i,t+1)		=	smap(etemp,utemp,dtemp);
		end
	end
end

peff						=	zeros(2,tg);
for t=1:tg
	peff(:,t)			=	pricespec002([0;0],[0;0],emap(smapinv(statecons(:,t),1))',dmap(dinit),umap(smapinv(statecons(i,t),2))',x,b,regrid,H,N,bmvargrid,0,0);
end


tvec	=	(1:tg-3);
f1		=	figure;
%  set(f1,'Renderer','painters');
%  set(f1, 'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 18 7.5]);




h1			=	line(tvec,pspec(1,3:tg-1),'Color','r');
ax1		=	gca;
title(['Positive Impulse repng:  \sigma_\eta^2 = ' num2str(pgrid(1),'%5.3f') '   \sigma_\gamma^2 = ' num2str(pgrid(2),'%5.3f') '  \rho_\eta  = ' num2str(pgrid(3),'%5.3f') '  \rho_\gamma = ' num2str(pgrid(4),'%5.3f') ' x = ' num2str(pgrid(5),'%5.3f') '  \sigma_x^2  = ' num2str(pgrid(6),'%5.3f')],'Interpreter','tex');
%  set(get(ax1,'Title'),'String',['Positive Impulse Simulation repng: \sigma_\eta^2 = ' num2str(pgrid(1),'%5.3f') ' | \sigma_\gamma^2 = ' num2str(pgrid(2),'%5.3f') ' | \rho_\eta = ' num2str(pgrid(3),'%5.3f') ' | \rho_\gamma = ' num2str(pgrid(4),'%5.3f') ' | \bar{x} ' num2str(pgrid(5),'%5.3f') ' | \sigma_x^2 = ' num2str(pgrid(6),'%5.3f')])
pos		=	get(ax1,'Position');

set(ax1,'XColor','k','YColor','k','Color','none','Position',pos);
set(get(ax1,'XLabel'),'String','Period');
set(get(ax1,'YLabel'),'String','Prices')




h2		=	line(tvec,pcons(1,3:tg-1),'Color','b','Parent',ax1);

h3		=	line(tvec,peff(1,3:tg-1),'Color','g','Parent',ax1);

h35	=	line(tvec,PriceSpecComp(1,3:tg-1),'Color','k','Parent',ax1);

ax2	=	axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','XAxisLocation','bottom','YColor','k','XColor','w');
h4		=	line(tvec,inv(1,3:tg-1),'Color',[0 1 1],'Parent',ax2);
h45	=	line(tvec,SpecInvComp(1,3:tg-1),'Color','m','Parent',ax2);
set(ax2,'YLim',[0 max(inv(1,:))*1.1+0.001], 'XLim',get(ax1,'XLim'))
set(get(ax2,'YLabel'),'String','Inventory')
legend([h1 h35 h2 h3 h4 h45],'Mon. Spec. Price', 'Comp. Spec. Price', 'Cons. Price', 'Eff. Price','Mon. Inventory','Comp. Inventory', 'Location','SouthEast');
print(f1,'-djpeg',['outgraph/posimpulse' int2str(repng) '.jpg']);
%  set(ax2,'Position',get(ax1,'Position'),'YAxisLocation','right','XAxisLocation','bottom');



%  set(get(ax3,'Title'),'String',['Negative Impulse Simulation repng: \sigma_\eta^2 = ' num2str(pgrid(1),'%5.4f') ' | \sigma_\gamma^2 = ' num2str(pgrid(2),'%5.4f') ' | \rho_\eta = ' num2str(pgrid(3),'%5.4f') ' | \rho_\gamma = ' num2str(pgrid(4),'%5.4f') ' | \bar{x} ' num2str(pgrid(5),'%5.3f') ' | \sigma_x^2 = ' num2str(pgrid(6),'%5.3f')])
f2			=	figure;
h5			=	line(tvec,pspec(2,3:tg-1),'Color','r');
ax3		=	gca;
title(['Negative Impulse repng:  \sigma_\eta^2  = ' num2str(pgrid(1),'%5.3f') '  \sigma_\gamma^2  = ' num2str(pgrid(2),'%5.3f') '  \rho_\eta  = ' num2str(pgrid(3),'%5.3f') '  \rho_\gamma  = ' num2str(pgrid(4),'%5.3f') ' x = ' num2str(pgrid(5),'%5.3f') ' \sigma_x^2 = ' num2str(pgrid(6),'%5.3f')],'Interpreter','tex');
pos		=	get(ax3,'Position');

set(ax3,'XColor','k','YColor','k','Color','none','Position',pos);
set(get(ax3,'XLabel'),'String','Period');
set(get(ax3,'YLabel'),'String','Prices')
h55	=	line(tvec,PriceSpecComp(2,3:tg-1),'Color','k','Parent',ax3);
hold all
h6		=	line(tvec,pcons(2,3:tg-1),'Color','b','Parent',ax3);
hold all
h7		=	line(tvec,peff(2,3:tg-1),'Color','g','Parent',ax3);
hold all
ax4	=	axes('Position',get(ax3,'Position'),'Color','none','YAxisLocation','right','XAxisLocation','bottom','YColor','k','XColor','w');

h8		=	line(tvec,inv(2,3:tg-1),'Parent',ax4,'Color',[0 1 1]);
h9		=	line(tvec,SpecInvComp(2,3:tg-1),'Color','m','Parent',ax4);
set(get(ax4,'YLabel'),'String','Inventory')
set(ax4,'YLim',[0 max(inv(2,:))*1.1+0.001], 'XLim',get(ax3,'XLim'))
legend([h5 h55 h6 h7 h8 h9],'Mon. Spec. Price', 'Comp. Spec. Price', 'Cons. Price', 'Eff. Price','Mon. Inventory', 'Comp. Inventory');

print(f2,'-djpeg',['outgraph/negimpulse' int2str(repng) '.jpg']);
save('outmat/impulserepng.mat','pspec','pcons','peff','PriceSpecComp','inv','SpecInvComp','tvec','tg','statespec','StateSpecComp','-V7.3');

quit
