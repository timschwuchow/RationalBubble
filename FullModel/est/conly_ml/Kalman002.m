%  Copyright 2012-2013 Timothy John Schwuchow
%  Kalman002.m		- Kalman Filter
%  Production Version

function [KalEst KalVar]		=	Kalman002(data,xiTTInit,PTTInit,H,F,Q)

T				=	numel(data);
xiTT_1		=	zeros(3,T);
xiTT_1(:,1)	=	xiTTInit;
PTT_1			=	cell(T,1);
PTT_1{1}		=	PTTInit;

KalEst		=	zeros(T,1);
KalVar		=	zeros(T,1);
KalEst(1)	=	data(1) - H'*xiTT_1(:,1);
KalVar(1)	=	H'*PTT_1{1}*H;

for t=1:T-1
	xiTT_1(:,t+1)		=	F*(xiTT_1(:,t) + PTT_1{t}*H*(H'*PTT_1{t}*H)^(-1)*(data(t)-H'*xiTT_1(:,t)));
	PTT_1{t+1}			=	F*(PTT_1{t} - PTT_1{t}*H*(H'*PTT_1{t}*H)^(-1)*H'*PTT_1{t})*F' + Q;
	KalEst(t+1)			=	(data(t+1)-H'*xiTT_1(:,t+1));
	KalVar(t+1)			=	H'*PTT_1{t+1}*H;
end
