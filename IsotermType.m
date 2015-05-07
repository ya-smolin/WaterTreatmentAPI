classdef IsotermType < handle
    
    properties
        name; %String
        lowerB; % double[]
        upperB; % double[]
        formula
    end
    
    methods
        %formula conventions
        %a = f(c)
        %k1,k2...k9 - koefitients
        %p1,p2...p9 - problem variables
        %use */+-^ don't use .*./.+.-.^
        function this = IsotermType(name, formula, lowerB, upperB)
            if nargin < 4
                this = 0;
                disp('incorrect number of arguments')
                return;
            end
            
            koefNames = unique(regexp(formula,'k\d', 'match')); %{p1,p2,..., p9}
            koefSizeMustBe = length(koefNames);
            koefSizeUserDef1 = length(lowerB);
            koefSizeUserDef2 = length(upperB);
            if (koefSizeUserDef1 ~= koefSizeMustBe || koefSizeUserDef2 ~= koefSizeMustBe)
                warning('bounds have incorect size correct size');
                lowerB = zeros(1, koefSizeMustBe);
                upperB = repmat(Isoterm.inf, 1, koefSizeMustBe);
            else
                %lowerB(lowerB == 0) = Isoterm.zero;
                %upperB(upperB == 0) = Isoterm.zero;
            end
            
            this.name = name;
            %lowerB(lowerB == 0) = Isoterm.zero;
            lowerB(lowerB == Inf) = Isoterm.inf;
            this.lowerB = lowerB;
            %upperB(upperB == 0) = Isoterm.zero;
            upperB(upperB == Inf) = Isoterm.inf;
            this.upperB = upperB;
            this.formula = formula;
        end
    end   
end

