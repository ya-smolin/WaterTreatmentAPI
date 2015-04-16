classdef IsotermModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
       zero = 1e-16;
    end
    properties
        isotermType; %IsotermType
        isotermResult;  %cfit
        gof; %struct( 'sse', [], 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    end
    
    methods
        function I = IsotermModel(isotermType, Cr, Ar)
            curIsotermFitmodel = isotermType.fitmodel;
            curConstants = isotermType.constants;
            lowerB = isotermType.lowerB;
            upperB = isotermType.upperB;
            
            %if bounds depends from Cr or Ar
            if(strcmp(isotermType.name,'Temkin'))
                lowerB = [0, 1/min(Cr)];
                %we do not have to use first point for temkin isoterm
            else
                Cr(end+1) = I.zero;
                Ar(end+1) = I.zero;        
            end
            
            isotermFun = isotermType.handler;
            pNum = length(probnames(curIsotermFitmodel));
            if ~isempty(curConstants)
                if(pNum == length(curConstants))
                    isotermFunPar = @(k)isotermFun(k, curConstants, Cr);
                else
                    error('pNum != length(curConstants)');
                end
            else
                if(pNum~=0)
                    disp('avto parameter fitting');
                    fstr = formula(isotermType.fitmodel);
                    kNum = numcoeffs(curIsotermFitmodel);
                    for j = 1:pNum
                        fstr = strrep(fstr, strcat('p',num2str(j)), strcat('k',num2str(kNum+1)));
                        lowerB(end+1) = max(Cr);
                        upperB(end+1) = Inf;
                    end
                    isotermType = IsotermType(isotermType.name, fstr, lowerB, upperB);
                    isotermFunPar = @(k)isotermType.handler(k, Cr);
                    curIsotermFitmodel = isotermType.fitmodel;
                else
                    isotermFunPar = @(k)isotermFun(k, Cr);
                end
            end
          
           
            errorIsotermFun = @(k)sum((isotermFunPar(k) - Ar).^2);
            
            opts = fitoptions(curIsotermFitmodel);
            opts.StartPoint = I.getStartPointUsingSA(errorIsotermFun, lowerB, upperB);
            opts.Display = 'Off';
            opts.Lower = lowerB;
            opts.Upper = upperB;
            opts.Robust = 'LAR';
            opts.Algorithm = 'Trust-Region';
            [xData, yData] = prepareCurveData(Cr, Ar);
            if ~isempty(curConstants)
                [fitResult, gof] = fit(xData,   yData,   curIsotermFitmodel, opts, 'problem', curConstants);
            else
                [fitResult, gof] = fit(xData,   yData,   curIsotermFitmodel, opts);
            end
            I.isotermType = isotermType;
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

