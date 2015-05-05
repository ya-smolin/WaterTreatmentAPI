classdef IsotermType < handle
    
    properties(Constant)
        zero = 1e-16;
    end
    
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
                lowerB = repmat(Isoterm.zero, 1, koefSizeMustBe);
                upperB = repmat(Isoterm.inf, 1, koefSizeMustBe);
            else
                lowerB(lowerB == 0) = this.zero;
                upperB(upperB == 0) = this.zero;
            end
            
            this.name = name;
            lowerB(lowerB == 0) = this.zero;
            this.lowerB = lowerB;
            upperB(upperB == 0) = this.zero;
            this.upperB = upperB;
            this.formula = formula;
        end
    end   
end

