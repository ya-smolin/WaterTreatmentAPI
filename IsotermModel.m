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
            isotermFun = isotermType.handler;
            curConstants = isotermType.constants;
            
            if(isempty(curConstants) && ~isempty(probnames(curIsotermFitmodel)))
                disp('you did not set isoterm constant parameter');
                return;
                %TODO: make avto replace parameter on variable
            end
            if ~isempty(curConstants)
                errorIsotermFun = @(k)sum((isotermFun(k, curConstants, Cr) - Ar).^2);  
            else
                errorIsotermFun = @(k)sum((isotermFun(k, Cr) - Ar).^2);
            end

            curLowerB = isotermType.lowerB;
            curUpperB = isotermType.upperB;

            opts = fitoptions(curIsotermFitmodel);
            populSize = 20;
            opts.StartPoint = I.getStartPointUsingGA(errorIsotermFun, numcoeffs(curIsotermFitmodel), curLowerB, curUpperB, populSize);
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
            options = gaoptimset(options,'PopInitRange', [0; 1000]);
            options = gaoptimset(options,'PopulationSize', populSize);
            options = gaoptimset(options,'Generations', 100);
            %options = gaoptimset(options,'StallGenLimit', 20);
            %options = gaoptimset(options,'TolFun', 0.0001);
            %options = gaoptimset(options,'Display', 'iter');
            %options = gaoptimset(options,'EliteCount', 3);
            %options = gaoptimset(options,'PlotFcns', { @gaplotbestf });
            startPoint = ga(fun, nvars,[],[],[],[],lb,ub,[],[],options);
        end
        
        function startPoint = getRandomStartPoint(dim, range)
            startPoint = random('Uniform', range(1), range(2), dim, 1);
        end
    end
    
end

