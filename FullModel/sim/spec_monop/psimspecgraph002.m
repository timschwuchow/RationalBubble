%  Copyright 2012 Timothy John Schwuchow
%  psimspecgraph002.m		-	Graph speculator and consumer-only models


clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');
gn		=	repng;
load(['outmat/pricesim' int2str(gn) '.mat']);


f1	=	figure;
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches') ;
set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
tvec	=	(1:tg);
h1	=	plot(tvec,pspec(1,1:tg),'-b');
hold on
h2	=	plot(tvec,pcons(1,1:tg),'-r');
hold on
h3	=	plot(tvec,peff(1,1:tg),'-g');
xlabel('Time');
ylabel('Prices');
legend('Speculator','Consumer','Efficient');
print(f1,'-djpeg',['../../graphs_t1/prices' int2str(gn) '.jpg']);
close(f1);

f2	=	figure;
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches') ;
set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
tvec	=	(1:tg);
h1	=	plot(tvec,inv(1,1:tg),'-b');
xlabel('Time');
ylabel('Inventory');
print(f2,'-djpeg',['../../graphs_t1/inv' int2str(gn) '.jpg']);
close(f2);

f3	=	figure;
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches') ;
set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
tvec	=	(1:tg);
h2	=	plot(tvec,umap(smapinv(statespec(1,1:tg),2)));
xlabel('Time');
ylabel('Common Innovation');
print(f3,'-djpeg',['../../graphs_t1/shock' int2str(gn) '.jpg']);
close(f3);

f4	=	figure;
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches') ;
set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
tvec	=	(1:tg);
h2	=	plot(tvec,dmap(smapinv(statespec(1,1:tg),3)), '-b');
hold on
h3	=	plot(tvec,dmap(smapinv(statecons(1,1:tg),3)), '-r');
xlabel('Time');
ylabel('Common Bias');
legend('Speculator', 'Consumer');
print(f4,'-djpeg',['../../graphs_t1/bias' int2str(gn) '.jpg']);
close(f4);

f5	=	figure;
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches') ;
set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
h2	=	plot(tvec,dmap(smapinv(statespec(1,1:tg),1)));
xlabel('Time');
ylabel('Common Flow Value');
print(f5,'-djpeg',['../../graphs_t1/eta' int2str(gn) '.jpg']);
close(f5);

f6	=	figure;
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches') ;
set(gcf, 'PaperPosition', [0.25, 1.5, 8.0, 8.0]);
h2	=	plot(tvec,cf(1,1:tg),'-r');
hold on
h3	=	plot(tvec,pi(1,1:tg),'-g');
xlabel('Time');
ylabel('Money');
legend('Discounted Cash Flows', 'Cumulative Profits');
print(f6,'-djpeg',['../../graphs_t1/cash' int2str(gn) '.jpg']);
close(f6);

quit