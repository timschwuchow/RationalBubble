%  Copyright 2012 Timothy John Schwuchow
%  psimspecdatdesc002.m		-	Describe features of LA data
clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');

indata 	= 	csvread('../../../data/flipmergeouttot001.csv');
price		=	indata(:,1);
period	=	indata(:,2);
inv		=	indata(:,3:4);

%  Plot price/quantity series
tvec		=	[1:numel(period)];
f1			=	figure;
h1			=	line(tvec,price,'Color','b');
ax1		=	gca;
set(get(ax1,'XLabel'),'String','Period');
set(get(ax1,'YLabel'),'String','Price Index');
set(ax1,'XColor','k','YColor','k','XLim',tvec([1 numel(period)]),'YLim',[min(price)*0.8 max(price)*1.2+0.0001]);
ax2 		=	axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right','Color','none','XColor','k','YColor','k','XLim',tvec([1 numel(period)]),'YLim',[min(inv(:,2))*0.8 max(inv(:,2))*1.2+0.0001]);
h2			=	line(tvec,inv(:,2),'Color','r','Parent',ax2);
set(get(ax2,'XLabel'),'String','Period');
set(get(ax2,'YLabel'),'String','Inventory');
legend([h1 h2],'Price Index', 'Inventory');
print(f1,'-djpeg','out1/laseries.jpg');
close(f1);
