classdef ZModel < handle
    properties
        d,
        al, be, Ar, Br, wr,
        Ed, A1, A2, A3, A4, B0, B1, St, D, D1, K1, K2,
        tau, nt, T, tDia, xDia, tInd,
        H, M, N, problemSize, split,
        q, Sf, S, Lf,
        Ss, Sfs, Sav,
        tGrid, rGrid, xGrid, uGrid
    end
    
    methods
        function o = ZModel()
            d = DataModule();
            o.d = d;
            %% Solution init size
            o.nt = 1; %observed cycles
            
            o.H = 9;  %x
            o.M = o.H; %r, %u
            o.N = 10;  %t
            
            o.al = 0;
            o.be = 2;
            
            o.tau = 1/(o.N-1);
            o.tGrid=0:o.tau:1;
            o.T = o.nt * d.tauRes;
            o.Ed = d.Ds * d.Dg * o.T / d.R^2;
            
            o.St = d.kf * o.T * (1 - d.eps) / (d.R * d.eps);
            o.D = d.Dg * d.Dx * o.T / d.L^2;
            o.B0 = d.Lmax / d.R;
            o.D1 = @(x, t) 3 * d.Dg * o.St * (1 + o.B0 * o.Lf.grid(x,t)).^2;
            
            o.A1 = d.K * o.tau * d.Xf * d.Lmax ./ (d.Y * d.q0 * d.ro * d.R);
            o.A2 = d.K * d.Xf * d.Lmax.^2 / (d.S0 * d.Df * d.Y);
            o.A3 = d.K * o.T;
            o.A4 = d.Kd * o.T;
            
            o.B1 = d.Ks / d.S0;
            o.K1 = o.T * d.Dg * d.Ap * d.kf * d.S0 / (12 * pi * d.Vc * d.ro * d.q0);
            o.K2 = o.T * d.Dg * d.K * d.Xf * d.Lmax / (4 * pi * d.R * d.ro * d.Y * d.q0);
            
            
            o.q = DFun([o.H, o.M, o.N]);
            o.S = DFun([o.H, o.N]);
            o.Sf = DFun([o.H, o.M, o.N]);
            o.Lf = DFun([o.H, o.N]);
            o.initConditions();
            
            o.Sfs = @(x,t)o.Sf.grid(x, o.M, t);
            o.Ss = @(x,t)o.Sf.grid(x, 1, t);
            o.Sav = @(x,t)sum(o.Sf.grid(x, :, t)) / (o.M * o.B1 + sum(o.Sf.grid(x, :, t)));
            
            %% Problem init size
            o.split = [o.q.sizeGridVar, o.S.sizeGridVar, o.Sf.sizeGridVar, o.Lf.sizeGridVar];
            o.problemSize = sum(o.split);
            o.tDia = o.q.varDia(o.q.dim);
            o.tInd = @(k)(k-1)*o.tau;
            o.xDia = o.q.varDia(1);
            
            %% Collocations init
            % create grid from roots of orthogonal polynomial
            rootsNum = o.M-2; %-1 for each boundary condition
            yP = YacobiPolynomial(rootsNum, o.al, o.be);
            o.rGrid = yP.rootsWithBoundary(1, 1);
            o.xGrid = o.rGrid;
            o.uGrid = o.rGrid;
            lPs = LagrangePolynomials(o.rGrid);
            % Ajm=Lm'(xj), Bjm=Lm''(xj)
            o.Ar = lPs.derValInNodes1(1);
            o.Br = lPs.derValInNodes1(2);
            o.wr = Radau.rdwExactInt(o.rGrid);
            % Unit tests section
            yP.unitTestOrthogonality();
            yP.unitTestYita();
            lPs.unitTestPolDer();
            lPs.unitTestAB();
            Radau.unitTestRDW();
            %% Solve model
            o.findSolution();
            o.test();
        end
        
        function initConditions(o)
            o.q.condInitial(0);
            o.S.condInitial(0);
            o.Lf.condInitial(o.d.L0 / o.d.Lmax);
            %o.Lf(:, o.N) = 1;
            o.S.condBoundLow(1, 1);
            %o.Sfs(1, :) = 1;
        end
        
        function X0 = getFirstApprox(o)
            X0 = zeros(o.problemSize, 1);
            splitInd = o.getSplitInd();
            X0(splitInd(4,1):splitInd(4,2)) = o.d.L0 / o.d.Lmax;
        end
        
        function Z = model(o, X, k)
            o.unpack(X, k);
            Z=[];
            for i = o.xDia
                Z(end+1) = eq_isoterm(o, i, k);
                for j = o.q.varDia(2)
                    Z(end+1) = eq_q(o, i, j, k);
                end
                for j = o.Sf.varDia(2)
                    Z(end+1) = eq_Sf(o, i, j, k);
                end
                Z(end+1) = eq_q_bound(o, i, k);
                Z(end+1) = eq_Lf(o, i, k);
                Z(end+1) = eq_S(o, i, k);
            end
        end
        
        function X0 = findSolution(o)
            options = optimoptions('fsolve', 'Display', 'iter',...
                'Algorithm', 'levenberg-marquardt', 'UseParallel', 1);
            %options.Algorithm = 'trust-region-reflective';
            X0 = getFirstApprox(o);
            for k = o.tDia
                X0 = fsolve(@(X)o.model(X, k), real(X0), options);
            end
        end
        
      function unpack(o, X, k)
        splitInd = o.getSplitInd();
        o.q.unpack(X(splitInd(1,1):splitInd(1,2)), k);
        o.S.unpack(X(splitInd(2,1):splitInd(2,2)), k);
        o.Sf.unpack(X(splitInd(3,1):splitInd(3,2)), k);
        o.Lf.unpack(X(splitInd(4,1):splitInd(4,2)), k);
      end
        
        
        function zero = eq_isoterm(o, i, k)
            zero = o.q.grid(i, o.M, k) - o.Ss(i, k).^o.d.n;
        end
        
        function zero = eq_S(o, i, k)
            if i == o.H
                zero = o.S.grid(i,k) - o.S.grid(i-1,k);
            else
                zero = -(1 / o.tau) * (o.S.grid(i,k) - o.S.grid(i,k-1)) + sum((o.D * o.Br(i,:)' - ...
                    o.nt * o.d.Dg * o.Ar(i,:)') .* o.S.grid(:,k)) - o.D1(i,k) * (o.S.grid(i,k) - o.Sfs(i,k));
            end
        end
        
        function zero = eq_q_bound(o, i, k)
            zero = (1 / o.tau) * dot(o.wr, o.q.grid(i, :, k) - o.q.grid(i, :, k-1)) - ...
                o.K1 * o.S.grid(i,k) + o.K1 * o.Sfs(i, k) + o.K2 * o.Lf.grid(i,k) * o.Sav(i, k);
        end
        
        function zero = eq_Lf(o, i, k)
            zero = o.Lf.grid(i, k) - o.Lf.grid(i, k-1) - o.tau * o.d.Dg * o.A3 * o.Sav(i, k) * ...
                o.Lf.grid(i,k) + o.tau * o.d.Dg * o.A4 * o.Lf.grid(i, k);
        end
        
        function zero = eq_Sf(o, i, j, k)
            zero = sum(o.Br(j, :) .* o.Sf.grid(i,:, k)) ...
                - o.A2 * o.Lf.grid(i, k)^2 * o.Sf.grid(i, j, k) / (o.B1 + o.Sf.grid(i, j, k));
        end
        
        function zero = eq_q(o, i, j, k)
            if j == 1
                %lower boundary
                zero = o.q.grid(i, 2, k) - o.q.grid(i, 1, k);
            else
                sumR = sum((o.d.s ./ o.rGrid(j) .* o.Ar(j, :) + o.Br(j, :)) .* o.q.grid(i, :, k));
                zero = o.q.grid(i, j, k-1) - o.q.grid(i, j, k) + o.tau * o.Ed * sumR;
            end
        end
        
        function splitInd = getSplitInd(o)
            funSize = length(o.split);
            splitInd = zeros(funSize, 2);
            for i = 1:funSize
                if i == 1
                    splitInd(1,1) = 1;
                else
                    splitInd(i,1) = splitInd(i-1, 2) + 1;
                end
                splitInd(i,2) = sum(o.split(1:i));
            end
        end
      
       function test(o)
            plotSize.M = 3;
            plotSize.N = 1;
            subplot(plotSize.M, plotSize.N, 1);
            axis([0, 1, 0, 1]);
            qSol = real(o.q.grid(o.H, :, o.N));
            q_val = reshape(qSol, 1, o.M);
            plot(o.rGrid, q_val);

            subplot(plotSize.M, plotSize.N, 2);
            axis([0, 1, 0, 1]);
            Lf_val = real(o.Lf.grid(o.H, :));
            plot(o.tGrid, Lf_val);

            subplot(plotSize.M,plotSize.N, 3);
            axis([0, 1, 0, 1]);
            SsLoc = reshape(real(o.Sf.grid(:, 1, :)), o.H, o.N);
            SfsLoc = reshape(real(o.Sf.grid(:, o.M, :)), o.H, o.N);
            plot(o.xGrid, SsLoc(:, o.N));
            hold on;
            plot(o.xGrid, SfsLoc(:, o.N));
        end
    end
end





