%  Copyright 2012 Timothy John Schwuchow
%  backestobj002.m	-	Test closeness of simulated correlation structure with true correlation structure

function obj	=	backestobj002(pg,b,T,H,N,x,lags,tdpcor)
[p eta delta inno]		=	pricesim002(pg,b,T,H,N,x);
[dpreg dpcor]				=	pricereggen002(p,lags);

cordif						=	dpcor - tdpcor;
obj							=	cordif*cordif';
