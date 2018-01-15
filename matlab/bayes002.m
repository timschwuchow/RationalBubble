%  Copyright 2012 Timothy John Schwuchow
%  bayes002.m	-	Solve for steady state variance of common/idiosyncratic bias and generate error coefficients used in model (correspond to c_{1\mu} terms and the like in theory write up

function [vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid mmgrid]		=	bayes002(pgrid,b,H,vz,q,T)
warning off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract simulation parameter vectors from grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vegrid		=	pgrid(:,1);
vggrid		=	pgrid(:,2);
regrid		=	pgrid(:,3);
rggrid		=	pgrid(:,4);
xvargrid	=	pgrid(:,6);
tolcrit		=	1e-8;
maxiter		=	10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create empty output vectors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n			=	numel(vegrid);
bmvargrid   =	zeros(size(vegrid));
vmgrid   =	zeros(size(vegrid));
vdgrid		=	zeros(size(vegrid));
vigrid		=	zeros(size(vegrid));
ce1grid		=	zeros(size(vegrid));
cg1grid		=	zeros(size(vegrid));
cd1grid		=	zeros(size(vegrid));
ci1grid		=	zeros(size(vegrid));
ce2grid		=	zeros(size(vegrid));
cg2grid		=	zeros(size(vegrid));
cd2grid		=	zeros(size(vegrid));
ci2grid		=	zeros(size(vegrid));
cb2grid		=	zeros(size(vegrid));
mmgrid		=	zeros(size(vegrid));

%  How this works: The outer loop cycles through each set of simulation parameters (from 1 to n) and makes initial guesses vd0 (the variance of delta_{t}) and vi0 (the variance of iota).
%  The inner loop finds the steady state variances of delta/iota by taking guesses vd0 and vi0 from the last iteration and computing the posterior variances (vd1, vi1) after one round of 'learning'.  The loop terminates when the difference between [vd0 vi0] and [vd1 vi1] is less than the tolerance.  While I cannot prove this loop always converges, in practice this hasn't been an issue.
for i=1:n
    i
	vd0		=	vegrid(i)/2;
	vi0		=	vggrid(i)/2;
    %guess at variance of marginal bid distribution
    vm0     =   vz/2; 
	tol		=	1;
	iter	=	0;
	while ( (tol > tolcrit) & (iter < maxiter) )
		iter					=	iter + 1;
		V						=	[vegrid(i) 0 0 0 0; 0 vggrid(i) 0 0 0; 0 0 vd0 0 0; 0 0 0 vi0 0; 0 0 0 0 vm0]; % Variance/covariance matrix of errors
		%  The last (empty) row of V is the 'variance' of the marginal bidder's idiosyncratic valuation (as perceived by consumers).  In both the consumer-only model and the naive-consumer with speculator model, consumers treat this as constant, hence the zero variance/covariance.  We'll need to modify this in the noise-trader model
		C1s						=	[1, 1, rggrid(i)-regrid(i), rggrid(i)-regrid(i),0]'; % Corresponds to c_1 vector in theory notes
		C1p						=	b*regrid(i)/(1-b*regrid(i))*[1, 0,-regrid(i), -regrid(i),0]'; % c_p in theory notes
		G1						=	V*C1s*(C1s'*V*C1s)^(-1)*C1s' - eye(size(V)); % Gamma_1 in theory notes
		clist1(1,:)				=	C1p'*G1; % Price coefficients (correspond to c_{1\mu} etc in equilibrium price equation from theory notes
		% Ex post signals
		C2s1					=	C1s;
		% Common-based signal
		C2s2					=	[1/(1-b*regrid(i)) + clist1(1,1),0,clist1(1,3) - regrid(i)/(1-b*regrid(i)),-regrid(i)/(1-b*regrid(i)),1]'; % c_2 in theory notes;;; why the -1 at the end?
		% Idiosyncratic-based signal
		%C2s3					=	[-1,clist1(1,2),regrid(i),regrid(i) + clist1(1,4),1]'; % Not in notes, but not needed
		C2s					=	[C2s2 C2s1];  % This corresponds to the stacked c_{12} matrix in the theory notes. Using any two vectors from the set (C2s1, C2s2, C2s3) delivers the same result.  Two vectors are used to capture the learning that takes place when consumers receive their valuation signals and again when prices are revealed - the choice of vectors is arbitrary because any two vectors are sufficient to describe what consumers 'learn' by the end of the period.
		G2						=	V*C2s*(C2s'*V*C2s)^(-1)*C2s' - eye(size(V)); % Gamma_2 in notes
		C2d					=	(1-b*regrid(i))*[-clist1(1,1),0,-clist1(1,3),0,-1]';  
		clist2(1,:)			=	C2d'*G2; % Coefficients that map realizations of errors today into delta/iota next period

		% Variances of delta_t/iota_t implied by coefficients, assuming the variance of delta_{t-1} and iota_{it-1} are vd0 and vi0
		vd1					= 	clist2(1,1)^2*vegrid(i) + clist2(1,3)^2*vd0+clist2(1,5)^2*vm0;
		vi1					=	clist2(1,2)^2*vggrid(i) + clist2(1,4)^2*vi0;
        
        %  Compute variance of bid distribution in consumer-only model
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
        bmvargrid(i)=	(xvargrid(i) + vggrid(i).*( (1+cg1grid(i)).^2+rggrid(i).^2./(1-rggrid(i).^2) +(ci1grid(i).*cg2grid(i)).^2./(1-ci2grid(i).^2) + 2*rggrid(i).*ci1grid(i).*cg2grid(i)./(1-rggrid(i).*ci2grid(i)))).^0.5;
        
        %calculate variance of marginal bid distribution
        %the mean and variance do not appear sensitive to the choice of
        %bmvargrid(i)
        [mm,vm1]=simvm(bmvargrid(i),H,vz,q);
        
        
        
		tol					=	max(abs([(vd1 - vd0) (vi1 - vi0) (vm1-vm0)]));
		vd0					=	vd1;
		vi0					=	vi1;
        vm0					=	vm1;
        
	end
	if tol < tolcrit
		% If convergence was successful, capture coefficients from last iteration
		vdgrid(i)		=	vd1;
		vigrid(i)		=	vi1;
        vmgrid(i)      =   vm1;
        mmgrid(i)      =   mm;
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


%  Compute variance of bid distribution in consumer-only model
bmvargrid				=	(xvargrid + vggrid.*( (1+cg1grid).^2+rggrid.^2./(1-rggrid.^2) +(ci1grid.*cg2grid).^2./(1-ci2grid.^2) + 2*rggrid.*ci1grid.*cg2grid./(1-rggrid.*ci2grid))).^0.5;

%  Bid distribution variance when consumers have perfect knowledge of state (so bid distribution variance is just the variance of idiosyncratic flow values)
varvgrid					=	(xvargrid + vggrid./(1-rggrid.^2)).^0.5;

