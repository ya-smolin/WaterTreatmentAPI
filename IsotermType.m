classdef IsotermType < handle
    
    properties(Constant)
        zero = 1e-16;
    end
    
    properties
        name; %String
        fitmodel; %fittype
        handler; %function handler
        lowerB; % double[]
        upperB; % double[]
        constants = [];
    end
    
    methods
        %formula conventions
        %a = f(c)
        %k1,k2...k9 - koefitients
        %p1,p2...p9 - problem variables
        %use */+-^ don't use .*./.+.-.^
        function I = IsotermType(name, formula, lowerB, upperB)
            if nargin < 4
                I = 0;
                disp('incorrect number of arguments')
                return;
            end
            
            params = unique(regexp(formula,'p\d', 'match'));
            if ~isempty(params)
              I.fitmodel = fittype(formula, 'independent', 'c', 'dependent', 'a', 'problem', params);
              formula = regexprep(formula,'(?<=p)\d','($0)');
              prefix = '@(k,p,c)';
            else
              I.fitmodel = fittype(formula, 'independent', 'c', 'dependent', 'a');  
              prefix = '@(k,c)';
            end

            %ki --> k(i)
            formula = regexprep(formula,'(?<=k)\d','($0)');
            %[*^\] --> [.*.^.\]
            formula = regexprep(formula,'[\^\*/]','.$0');
            func = strcat(prefix, formula);
            I.handler = str2func(func);
            
            I.name = name;
            lowerB(lowerB == 0) = I.zero;
            I.lowerB = lowerB;
            upperB(upperB == 0) = I.zero;
            I.upperB = upperB;
        end
        
        function set.constants(I, constants)
            numPar = length(probnames(I.fitmodel));
            if numPar == 0
                disp('there are no any problem parameters, you should not use this method at all. Be calm and ride bike');
            end
            if length(constants) ~= numPar
                disp('expected ', numPar, ' but you seted ', len(constants))
            else
                %here code place for adequate people
                I.constants = constants;
            end
        end
        
        function constants = get.constants(I)
            constants = I.constants;
        end
    end   
end

