%  Copyright 2012-2013 Timothy John Schwuchow
%  bayes002_WV.m	-	Test build of Bayesian learning protocol with inventory constant (reparameterization consistent with brstheory_eq001.pdf theory results)
%  Working Version
function [vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid]		=	bayes002_WV(pgrid,b)

vegrid	=	pgrid(:,1);
vggrid	=	pgrid(:,2);
regrid	=	pgrid(:,3);
rggrid	=	pgrid(:,4);
xvargrid	=	pgrid(:,6);
tolcrit	=	1e-14;
maxiter	=	15000;


n			=	numel(vegrid);
vdgrid	=	zeros(size(vegrid));
vigrid	=	zeros(size(vegrid));
ce1grid	=	zeros(size(vegrid));
cg1grid	=	zeros(size(vegrid));
cd1grid	=	zeros(size(vegrid));
ci1grid	=	zeros(size(vegrid));
ce2grid	=	zeros(size(vegrid));
cg2grid	=	zeros(size(vegrid));
cd2grid	=	zeros(size(vegrid));
ci2grid	=	zeros(size(vegrid));
cb2grid	=	zeros(size(vegrid));
clist2	=	ones(1,5);
clist1	=	ones(1,5);
warning off
for i=1:n
	vd0		=	vegrid(i)/2;
	vi0		=	vggrid(i)/2;
	tol1		=	1;
	tol2		=	2;
	iter		=	0;
	while ( (tol2 > tolcrit | tol1 > tolcrit) & (iter < maxiter) )

		clist1old			=	clist1;
		clist2old			=	clist2;
		iter					=	iter + 1;
		V						=	[vegrid(i) 0 0 0 0; 0 vggrid(i) 0 0 0; 0 0 vd0 0 0; 0 0 0 vi0 0; 0 0 0 0 0];
		CSig1					=	[1, 1, rggrid(i)-regrid(i), rggrid(i)-regrid(i),0]';
		CExpBid				=	[b*clist1(1,3)*clist2(1,1) + b*regrid(i)/(1-b*regrid(i)),0,b*clist1(1,3)*clist2(1,3) - b*regrid(i)^2/(1-b*regrid(i)),- b*regrid(i)^2/(1-b*regrid(i)),0]';

		CFixBid				=	[-b*regrid(i)/(1-b*regrid(i)),0,b*regrid(i)^2/(1-b*regrid(i)),b*regrid(i)^2/(1-b*regrid(i)),0]';

%  		C1p					=	b*regrid(i)/(1-b*regrid(i))*[1, 0,-regrid(i), -regrid(i),0]';
%  		G1						=	V*CSig1*(CSig1'*V*CSig1)^(-1)*CSig1' - eye(size(V));
%  		clist1(1,:)			=	C1p'*G1;
%  		clistexp1(1,:)

		clist1(1,:)		=	CExpBid'*V*CSig1*(CSig1'*V*CSig1)^(-1)*CSig1' + CFixBid';
		% Ex post signals
		C2s1					=	CSig1;
		% Common-based signal
%  		CSig2					=	[-1/(1-b*regrid(i)) - clist1(1,1),0,-clist1(1,3) + regrid(i)/(1-b*regrid(i)),+regrid(i)/(1-b*regrid(i)),-1]';
		CSig2					=	[1/(1-b*regrid(i)) + clist1(1,1),0,clist1(1,3) - regrid(i)/(1-b*regrid(i)),-regrid(i)/(1-b*regrid(i)),-1]';
		% Idiosyncratic-based signal
		CSig3					=	[-1,clist1(1,2),regrid(i),regrid(i) + clist1(1,4),1]';
		C2s					=	[CSig2 CSig3];
		G2						=	V*C2s*(C2s'*V*C2s)^(-1)*C2s' - [1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0];
%  		G2						=	C2s*C2s'
		C2d					=	(1-b*regrid(i))*[-clist1(1,1),0,-clist1(1,3),0,-1]';
		clist2(1,:)			=	C2d'*G2;

%  		[clist1;clist2]
		vd1					= 	clist2(1,1)^2*vegrid(i) + clist2(1,3)^2*vd0;
		vi1					=	clist2(1,2)^2*vggrid(i) + clist2(1,4)^2*vi0;
		tol1					=	max(abs([(vd1 - vd0) (vi1 - vi0)]));
		tol2					=	max(abs([clist1 clist2] - [clist1old clist2old]));
%  		fprintf('Iter %d | tol %5.4g\n',iter,tol);
		vd0					=	vd1;
		vi0					=	vi1;
%  		[tol vd1 vi1]
%  		fprintf('Iter %d Tol1 %5.4g Tol2 %5.4g \nCoefficients 1: %5.3f %5.3f %5.3f %5.3f %5.3f\nCoefficients 2: %5.3f %5.3f %5.3f %5.3f %5.3f\n******************\n',iter,tol1,tol2,clist1,clist2)
%  		fprintf('Tol1 %5.3g Tol2 %5.3g Iter %d\n',tol1,tol2,iter)
	end
%  	iter
%  	if mod(iter,5)==0
%  	end
	if ( (tol1 < tolcrit) & (tol2 < tolcrit) )
		vdgrid(i)		=	vd1;
		vigrid(i)		=	vi1;
		ce1grid(i)		=	clist1(1,1);
		cg1grid(i)		=	clist1(1,2);
		cd1grid(i)		=	clist1(1,3);
		ci1grid(i)		=	clist1(1,4);
		ce2grid(i)		=	clist2(1,1);
		cg2grid(i)		=	clist2(1,2);
		cd2grid(i)		=	clist2(1,3);
		ci2grid(i)		=	clist2(1,4);
		cb2grid(i)		=	clist2(1,5);
	else
		fprintf('No convergence for ve = %7.5f | vg = %7.5f | re = %7.5f | rg = %7.5f \n iter = %d Tol1 = %5.4g Tol2 = %5.4g | vdlast = %7.5f \n',vegrid(i),vggrid(i),regrid(i),rggrid(i),iter,tol1,tol2,vd1);
	end
end


bmvargrid				=	(xvargrid + vggrid.*( (1+cg1grid).^2+rggrid.^2./(1-rggrid.^2) +(ci1grid.*cg2grid).^2./(1-ci2grid.^2) + 2*rggrid.*ci1grid.*cg2grid./(1-rggrid.*ci2grid))).^0.5;

varvgrid					=	(xvargrid + vggrid./(1-rggrid.^2)).^0.5;

