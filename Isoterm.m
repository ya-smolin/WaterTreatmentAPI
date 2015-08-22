classdef Isoterm < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        zero = 1e-9;
        inf = 1e7;
    end
    
    properties
        id;
        isotermType;
        isotermResult;  %cfit
        gof; %struct( 'sse', [], 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    end
    
    methods
        %formula conventions
        %a = f(c)
        %k1,k2...k9 - koefitients
        %p1,p2...p9 - problem variables
        %use */+-^ don't use .*./.+.-.^
        function this = Isoterm(isotermType, Cr, Ar, constants)
            %==check input
            if nargin < 3
                error('uncorrect input argummets number')
            elseif nargin == 3
                constants = [];
            end
            
            name = isotermType.name;
            formula = isotermType.formula;
            lowerB = isotermType.lowerB;
            upperB = isotermType.upperB;
            
            %=resolve koef issue
            koefNames = unique(regexp(formula,'k\d', 'match')); %{p1,p2,..., p9}
            koefSizeMustBe = length(koefNames);
            %=resolve constants issue
            constNames = unique(regexp(formula,'p\d', 'match')); %{p1,p2,..., p9}
            constantSizeMustBe = length(constNames);
            contantSizeUserDef = length(constants);
            if constantSizeMustBe == 0
                if contantSizeUserDef > 0
                    warning('constantSizeMustBe == 0, but contantSizeUserDef > 0');
                end
                %everything is OK. No constants
            else
                if constantSizeMustBe == contantSizeUserDef
                    %everything is OK. Contants defined correctly
                else
                    warning('autofitting all parameters\n constantSizeMustBe ~= contantSizeUserDef');
                    %autofit all parameters
                    %make all constants parameters
                    for j = 1:constantSizeMustBe
                        constants = [];
                        constNames = {};
                        formula = strrep(formula, strcat('p',num2str(j)), strcat('k',num2str(koefSizeMustBe+1)));
                        lowerB(end+1) = 0;
                        upperB(end+1) = Isoterm.inf;
                    end
                end
            end
            
            %if bounds depends from Cr or Ar
            if(strcmp(name,'Temkin'))
                lowerB = [0, 1/min(Cr)];
                %we do not have to use first point for temkin isoterm
            else
                Cr(end+1) = Isoterm.zero;
                Ar(end+1) = Isoterm.zero;
            end
            
            fitmodel = fittype(formula, 'independent', 'c', 'dependent', 'a', 'problem', constNames);
            hFun = Isoterm.strFormulaToHandle(formula, constNames);
            isotermFunPar = @(k)hFun(k, Cr);
            errorIsotermFun = @(k)sum((isotermFunPar(k) - Ar).^2);
            opts = fitoptions(fitmodel);
%             opts.StartPoint = Isoterm.getStartPointUsingSA(errorIsotermFun, lowerB, upperB);
            opts.StartPoint = Isoterm.getStartPointUsingUniform(lowerB, upperB);
            opts.Display = 'Off';
            opts.Lower = lowerB;
            opts.Upper = upperB;
            opts.Robust = 'LAR';
            opts.Algorithm = 'Trust-Region';
            [xData, yData] = prepareCurveData(Cr, Ar);
            if ~isempty(constants) && ~isempty(constNames)
                [isotermResult, gof] = fit(xData,   yData,   fitmodel, opts, 'problem', constants);
            else
                [isotermResult, gof] = fit(xData,   yData,   fitmodel, opts);
            end
            
            this.isotermType = IsotermType(name, formula, lowerB, upperB);
            this.isotermResult = isotermResult;
            this.gof = gof;
        end
    end
    
    methods(Static)
        function hFun = strFormulaToHandle(formula, params)
            if ~isempty(params)
                formula = regexprep(formula,'(?<=p)\d','($0)');
                prefix = '@(k,p,c)';
            else
                prefix = '@(k,c)';
            end
            %ki --> k(i)
            formula = regexprep(formula,'(?<=k)\d','($0)');
            %[*^\] --> [.*.^.\]
            formula = regexprep(formula,'[\^\*/]','.$0');
            formula = strcat(prefix, formula);
            hFun = str2func(formula);
        end  
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
            %disp(opts.StartPoint)
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
            x0 = Isoterm.getStartPointUsingUniform(lb, ub);
            
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

