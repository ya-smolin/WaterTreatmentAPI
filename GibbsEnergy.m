classdef GibbsEnergy
    %GIBBSENERGY for adsortion in water
    %INPUT
    %   Cr: לד/כ
    %   Ar: לד/כ
    %   Va: סל^3/ד
    %   T:  Ê
    %   M:  לד/ללמכü
    %   Vz: סל^3/ללמכü
    %OUTPUT
    %   dG: ÊÄז/למכü
    %   lnK: log of ctg of line slope
    
    
    
    properties(Constant)
        MILLI = 1e3;
        R = 8.31/GibbsEnergy.MILLI;
        C_H2O = 55555;
        V_H2O = 0.018;
    end
    
    methods(Static)
        %   fitResult: x: 10^3/C_i דÎ/כ
        %              y: Va/a^*_i סל^3/לדÎ
        function [dG, lnK, fitResult, gof] = calculate_xpk(Cr, Ar, Va, T)
            x=GibbsEnergy.MILLI/Cr; %x
            y=Va./Ar; %y
            [fitResult, gof]=fit(x,y,'poly1');
            koef = coeffvalues(fitResult);
            k = koef(1);
            ctg=1./k;
            lnK=log(ctg);
            dG=GibbsEnergy.R*T*lnK;
        end
        
        %   fitResult: x: lg(y/x) 
        %              y: \Theta
        function [dG, lnK, fitResult, gof] = calculate(Cr, Ar, Va, T, Vz, M)
            Cr = Cr ./ M;
            Ar = Ar ./ M;%ללמכü/ךד
            x=Cr./(Cr + GibbsEnergy.C_H2O);
            A_H2O=(Va-Ar*Vz)/GibbsEnergy.V_H2O;
            y=Ar./(Ar+A_H2O);
            Teta=Ar*Vz/Va; %x
            Y=log10(y./x); %y
            [fitResult, gof]=fit(Teta,Y,'poly1');
            lgK=fitResult(0);
            K=10^lgK;
            lnK=log(K);
            dG=GibbsEnergy.R*T*lnK;
        end
        
        
    end
    
    
end

