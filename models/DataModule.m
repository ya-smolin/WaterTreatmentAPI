classdef DataModule < handle
    properties
        KF, n, s, Ds, kf, K, Ks, Y, Kd, Xf, L0, Lmax, S0, L, ro, R, Dx, Df, 
        eps, Vx, tauRes, q0, Dg, Ap, Vc
    end
    
    methods
        function o = DataModule()
            %% expiremental data in [gr cm min]
            MGL2GRCM3 = 10^(-3)/10^3;
            o.KF=11.48;
            o.n = 0.423;
            o.s = 2;
            o.Ds=1e-8 *60;%cm^2/s-->cm^2/min
            o.kf=6.7e-4 *60; %cm/s-->cm/min
            o.K =6.67e-3;%1/min
            o.Ks=125 * MGL2GRCM3; %mg/l-->gr/cm^3
            o.Y=1.35;%mg/mg-->gr/gr
            o.Kd=5e-5; %1/min
            o.Xf = 7.6* 10^(-3); %mg/cm^3-->gr/cm^3
            
            o.L0 = 5e-4;%cm
            o.Lmax= 1e-2;%cm
            o.S0=81* MGL2GRCM3; %mg/l-->gr/cm^3
            o.L=61.75; %XL cm
            o.ro=0.67; %gr/cm^3
            o.R=2.75e-2;%RAD cm

            W=350; %gr
            DIA = 5.08; %cm
            Q=234; %Q ml/min-->cm^3/min
           
            o.Dx=0 *60;%Dd cm^2/s-->cm^2/min
            o.Df=4.48e-6 *60;%cm^2/s-->cm^2/min
            %% derived data
            Vc=2*pi*DIA*o.L;
            roReal = W/Vc;
            o.eps=(o.ro - roReal)/o.ro;
            A = pi * (DIA/2)^2;
            o.Vx=Q/(o.eps*A);
            o.tauRes = o.L / o.Vx;
            o.q0=o.KF*o.S0^o.n; %my guess
            o.Dg=o.ro*o.q0*(1-o.eps)/(o.eps*o.S0);
            
            o.Ap = 4*pi*o.R^2; %hardcode for s=2
            o.Vc = 4/3*pi*o.R^3; %hardcode for s=2;
        end
        
    end
end

