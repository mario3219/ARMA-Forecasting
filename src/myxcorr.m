function myxcorr(y,x,noLags)

figure
[Cxy,lags] = xcorr( y, x, noLags, "coeff");
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(y) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off