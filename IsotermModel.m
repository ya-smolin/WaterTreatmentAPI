classdef IsotermModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        isotermResult;  %cfit
        gof = struct( 'sse', [], 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    end
    
    methods
        function I = IsotermModel(isotermType, Cr, Ar)
            curIsotermFitmodel = isotermType.fitmodel;
            curConstants = isotermType.constants;
            
            if(isempty(curConstants) && ~isempty(probnames(curIsotermFitmodel)))
                disp('you did not set isoterm constant parameter');
                return;
            end
            
            isotermFun = isotermType.handler;
            if ~isempty(curConstants)
                errorIsotermFun = @(k)sum((isotermFun(k, curConstants, Cr) - Ar).^2);  
            else
                errorIsotermFun = @(k)sum((isotermFun(k, Cr) - Ar).^2);
            end

            curLowerB = isotermType.lowerB;
            curUpperB = isotermType.upperB;

            opts = fitoptions(curIsotermFitmodel);
            populSize = 4 ^ length(curUpperB);
            %opts.StartPoint = I.getStartPointUsingGA(errorIsotermFun, numcoeffs(curIsotermFitmodel), curLowerB, curUpperB, populSize);
            opts.StartPoint = I.getStartPointUsingSA(errorIsotermFun, curLowerB, curUpperB);
            opts.Display = 'Off';
            opts.Lower = curLowerB;
            opts.Upper = curUpperB;
            opts.Robust = 'LAR';
            opts.Algorithm = 'Trust-Region';
            [xData, yData] = prepareCurveData(Cr, Ar);
            if ~isempty(curConstants)
                [fitResult, gof] = fit(xData,   yData,   curIsotermFitmodel, opts, 'problem',curConstants);
            else
                [fitResult, gof] = fit(xData,   yData,   curIsotermFitmodel, opts);
            end
            I.isotermResult = fitResult;
            I.gof = gof;
        end    
    end
    
      methods(Static)
        function startPoint = getStartPointUsingGA(fun, nvars, lb, ub, populSize)
            %% Start with the default options
            options = gaoptimset;
            %% Modify options setting
            %options = gaoptimset(options,'FitnessLimit', 0.8);
            %options = gaoptimset(options,'InitialPopulation', initialPopulation);
            options = gaoptimset(options,'PopInitRange',  [0; 100]);
            options = gaoptimset(options,'PopulationSize', populSize);
            options = gaoptimset(options,'Generations', 100);
            %options = gaoptimset(options,'StallGenLimit', 20);
            %options = gaoptimset(options,'TolFun', 0.0001);
            %options = gaoptimset(options,'Display', 'iter');
            %options = gaoptimset(options,'EliteCount', 3);
            %options = gaoptimset(options,'PlotFcns', { @gaplotbestf });
            startPoint = ga(fun, nvars,[],[],[],[],lb,ub,[],[],options);
        end
        function x0 = getStartPointUsingUniform(lb, ub)
            x0 = random('Uniform', lb, ub);
        end
        
        function startPoint = getStartPointUsingSA(fun, lb, ub)
            K = length(lb);
            N = length(ub);
            if K~=N
                error('the length of upper bound does not equal lower.');
            end
            %% This is an auto generated MATLAB file from Optimization Tool.
            x0 = IsotermModel.getStartPointUsingUniform(lb, ub);
            
            %% Start with the default options
            options = saoptimset;
            %% Modify options setting
            options = saoptimset(options,'Display', 'off');
            options = saoptimset(options,'HybridInterval', 'end');
            startPoint = simulannealbnd(fun,x0,lb,ub,options);
            display(startPoint);
        end

      end
      
end

