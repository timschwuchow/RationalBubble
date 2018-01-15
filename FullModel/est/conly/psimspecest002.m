%  Copyright 2012 Timothy John Schwuchow
%  psimspecest002.m		-	Estimate speculation model using LA data
clear all
format short g
path(path,'../../functions');
path(path,'../../functions/tests');


%  load('../../sim/t4/outmat/pricesim1.mat');

indata 	= 	csvread('../../../data/flipmergeouttot001.csv');
price		=	indata(:,1);
period	=	indata(:,2);
inv		=	indata(:,3:4);

%  Simulation parameters
vegss		=	2.00;
vggss		=	5.0;
regss		=	0.90;
rggss		=	0.2;
pargss	=	[vegss vggss regss rggss];
parub		=	[40.0 40.0 0.995 0.995];
parlb		=	[0.01 0.01 0.001 0.001];


%  Other parameters
x		=	1.0;
xvar	=	2.5;
H		=	0.45;
N		=	1;
b		=	0.99;

%  Options
optcon	=	optimset('Display', 'iter','TolFun',1e-8,'Algorithm','interior-point');
opt	=	optimset('Display', 'iter','TolFun',1e-8,'TolX',1e-12);
[parest1 obj1]	=	fmincon( @(par) psimspecestobj002(par,price,inv(:,2),zeros(size(inv(:,2))),x,xvar,H,N,b),pargss,[],[],[],[],parlb,parub,[],optcon);
[parest2 obj2]	=	fminsearch( @(par) psimspecestobj002(par,price,inv(:,2),zeros(size(inv(:,2))),x,xvar,H,N,b),pargss,opt);

[parest1 obj1; parest2 obj2]

