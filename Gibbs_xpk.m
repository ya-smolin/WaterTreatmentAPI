function [dG, lnK]=Gibbs_xpk(Cr,Ar,Va,T)
ro = 1e3;
R=8.31;

x=ro./Cr;
y=Va./Ar;
[fitResult, gof]=fit(x,y,'poly1');
koef = coeffvalues(fitResult);
k = koef(1);
ctg=1./k;
lnK=log(ctg); %K=ctg/
dG=(R/ro)*T*lnK;

xlabel('10^3/C_i ד־/כ');
ylabel('Va/a^*_i סל^3/לד־');

plot(x,y,'ko');
hold on;
plot(fitResult);
    leg=legend('experimantal points(10^3/C_i ; Va/a^*_i) ',...
    [' Gibbs free energy \DeltaG=',num2str(dG),' Êִז/למכü',...
    ' R^2=',num2str(gof.rsquare),'  lnK=', num2str(lnK)]);
    set(leg, 'Location','SouthEast');   
end

