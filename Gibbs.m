%%Gibbs free energy  H2O 20^oC
function [dG, lnK]=Gibbs(Cr,Ar,M,Va,Vz,T)
C_H2O=55555;
R=8.31;
V_H2O=0.018;

Cr=Cr./M;
Ar=Ar./M;%ללמכü/ךד
x=Cr./(Cr+C_H2O);
A_H2O=(Va-Ar*Vz)/V_H2O;
y=Ar./(Ar+A_H2O);
Y=log10(y./x);

Teta=Ar*Vz/Va;


[fitResult, gof]=fit(Teta,Y,'poly1');

lgK=fitResult(0);
K=10^lgK;
lnK=log(K);

dG=R/1000*T*lnK;

plot(Teta,Y,'ko');
hold on;
fplot(fitResult, [0, max(Teta)]);
leg=legend('experimantal points[lg(y/x) ; \Theta]',...
    [' Gibbs free energy \DeltaG=',num2str(dG),' ÊÄז/למכü',...
    ' R^2=',num2str(gof.rsquare),'  lnK=', num2str(lnK)]);
xlabel('lg(y/x) ללמכü/ד');
ylabel('\Theta ללמכü/ד');
set(leg, 'Location','SouthEast');
end