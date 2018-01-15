%  Copyright 2012-2013 Timothy John Schwuchow
%  Kalman002_WV.m		- Kalman Filter
%  Working Version

function [KalEst KalVar]		=	Kalman002_WV(data,XiTTInitDraws,PTTInit,Psi,F,Q)

T					=	numel(data);
ndraws			=	size(XiTTInitDraws,2);
CellXiTT			=	cell(T,1);


CellXiTT{1}		=	XiTTInitDraws;

PTT_1				=	cell(T,1);
PTT_1{1}			=	PTTInit;
KalEst			=	zeros(T,ndraws);
KalVar			=	zeros(T,1);
KalEst(1,:)		=	data(1) - Psi'*CellXiTT{1};
KalVar(1)		=	Psi'*PTT_1{1}*Psi;
for t=1:T-1
	CellXiTT{t+1}		=	F*(CellXiTT{t} + PTT_1{t}*Psi*(KalVar(t))^(-1)*(KalEst(t,:)));
	PTT_1{t+1}			=	F*(PTT_1{t} - PTT_1{t}*Psi*(KalVar(t))^(-1)*Psi'*PTT_1{t})*F' + Q;
	KalEst(t+1,:)		=	(data(t+1)-Psi'*CellXiTT{t+1});
	KalVar(t+1)			=	Psi'*PTT_1{t+1}*Psi;
end
