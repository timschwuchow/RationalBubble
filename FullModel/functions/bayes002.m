%  Copyright 2012 Timothy John Schwuchow
%  bayes002.m	-	Test build of Bayesian learning protocol with inventory constant (reparameterization consistent with brstheory_eq001.pdf theory results)

function [vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid]		=	bayes002(pgrid,b)

vegrid	=	pgrid(:,1);
vggrid	=	pgrid(:,2);
regrid	=	pgrid(:,3);
rggrid	=	pgrid(:,4);
xvargrid	=	pgrid(:,6);
tolcrit	=	1e-10;
maxiter	=	10000;


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
warning off
for i=1:n
	vd0		=	vegrid(i)/2;
	vi0		=	vggrid(i)/2;
	tol		=	1;
	iter		=	0;
	while ( (tol > tolcrit) & (iter < maxiter) )
		iter					=	iter + 1;
		V						=	[vegrid(i) 0 0 0 0; 0 vggrid(i) 0 0 0; 0 0 vd0 0 0; 0 0 0 vi0 0; 0 0 0 0 0];
		C1s					=	[1, 1, rggrid(i)-regrid(i), rggrid(i)-regrid(i),0]';
		C1p					=	b*regrid(i)/(1-b*regrid(i))*[1, 0,-regrid(i), -regrid(i),0]';
		G1						=	V*C1s*(C1s'*V*C1s)^(-1)*C1s' - eye(size(V));
		clist1(1,:)			=	C1p'*G1;
		% Ex post signals
		C2s1					=	C1s;
		% Common-based signal
		C2s2					=	[1/(1-b*regrid(i)) + clist1(1,1),0,clist1(1,3) - regrid(i)/(1-b*regrid(i)),-regrid(i)/(1-b*regrid(i)),-1]';
		% Idiosyncratic-based signal
		C2s3					=	[-1,clist1(1,2),regrid(i),regrid(i) + clist1(1,4),1]';
		C2s					=	[C2s2 C2s3];
		G2						=	V*C2s*(C2s'*V*C2s)^(-1)*C2s' - [1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0];
%  		G2						=	C2s*C2s'
		C2d					=	(1-b*regrid(i))*[-clist1(1,1),0,-clist1(1,3),0,-1]';
		clist2(1,:)			=	C2d'*G2;

%  		[clist1;clist2]
		vd1					= 	clist2(1,1)^2*vegrid(i) + clist2(1,3)^2*vd0;
		vi1					=	clist2(1,2)^2*vggrid(i) + clist2(1,4)^2*vi0;
		tol					=	max(abs([(vd1 - vd0) (vi1 - vi0)]));
%  		fprintf('Iter %d | tol %5.4g\n',iter,tol);
		vd0					=	vd1;
		vi0					=	vi1;
%  		[tol vd1 vi1]
	end
	if tol < tolcrit
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
		fprintf('No convergence for ve = %7.5f | vg = %7.5f | re = %7.5f | rg = %7.5f \n Tol = %8.5g | vdlast = %7.5f \n',vegrid(i),vggrid(i),regrid(i),rggrid(i),tol,vd1);
	end
end



bmvargrid				=	(xvargrid + vggrid.*( (1+cg1grid).^2+rggrid.^2./(1-rggrid.^2) +(ci1grid.*cg2grid).^2./(1-ci2grid.^2) + 2*rggrid.*ci1grid.*cg2grid./(1-rggrid.*ci2grid))).^0.5;

varvgrid					=	(xvargrid + vggrid./(1-rggrid.^2)).^0.5;

