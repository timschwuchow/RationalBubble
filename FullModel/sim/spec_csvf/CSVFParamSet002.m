%  Copyright 2012-2013 Timothy John Schwuchow
%  CSVFParamSet002.m		-	Set CSVF Parameters
%  Production  version
function [cd1grid,ce1grid,cd2grid,ce2grid,cb2grid,bmvargrid,imap,emap,dmap,umap,smap,smapinv] = CSVFParamSet002(ve,vg,re,rg,x,xvar,b,vdiag,H,N,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,ct)
path(path,'../../functions');
%  Grid out parameter combinations
[vegrid vggrid regrid rggrid xgrid xvargrid pgrid ngrid]		=	grid002(ve,vg,re,rg,x,xvar,b,vdiag);
clear ve vg re rg x xvar
%  Generate constants
[vdgrid vigrid ce1grid cg1grid cd1grid ci1grid ce2grid cg2grid cd2grid ci2grid cb2grid bmvargrid varvgrid]		=	bayes002_WV(pgrid,b);



%  Construct state grid
[imap emap umap dmap smap smapinv nstate]	=	statemappar002(pgrid,vdgrid,cb2grid,cd2grid,bmvargrid,ipts,epts,upts,dpts,maxshare,esd,usd,dsd,H,N);

if not(exist('outcsv/parameters.csv')==2)
	outfile = fopen('outcsv/parameters.csv','w');
	fprintf(outfile,'sim \tve \tvg \tre \trg \tx \tvx \tH \tb \tepts \tupts \tdpts \tipts \tvd \tvi \tce1 \tcd1 \tce2 \tcd2 \tcb2 \tbmvar \n');
	fclose(outfile);
end
outfile = fopen('outcsv/parameters.csv','a');
for i=1:ngrid
	fprintf(outfile,['%d' repmat('\t%5.3g',[1 size([pgrid  H b epts upts dpts ipts vdgrid vigrid ce1grid cd1grid ce2grid cd2grid cb2grid bmvargrid],2)]) '\n'],[ct pgrid(i,:) H b epts upts dpts ipts vdgrid(i) vigrid(i) ce1grid(i) cd1grid(i) ce2grid(i) cd2grid(i) cb2grid(i) bmvargrid(i)]);
end
fclose(outfile);
