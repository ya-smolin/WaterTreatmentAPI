classdef IsotermTableRow
    %TABLEROWVIEW Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        columnName =      {'name',    'formula',  'koef',     'params',   'show',     'rsquare', 'sse',     'koef range'};
        columnFormat =    {'char',    'char',     'char',  'char',  'logical',  'numeric', 'numeric', 'char'}
        columnEditable =  [false,   false,    false,    true,     true,     false,   false,   false]
        size = length(IsotermTableRow.columnName);
        CONFIDENCE_LEVEL_DEFAULT = 0.95;
    end
    properties
        data;
    end
    
    methods
        function this = IsotermTableRow(isoterm, confidenceLevel)
            if nargin == 1
                confidenceLevel = IsotermTableRow.CONFIDENCE_LEVEL_DEFAULT;
            end
            if isa(isoterm, 'IsotermType')
                this.data = {
                    isoterm.name,...
                    isoterm.formula,...
                    '',...
                    '',...
                    false,...
                    0,...
                    Inf,...
                    ''
                    };
            else
                isotermType = isoterm.isotermType;
                isotermResult =  isoterm.isotermResult;
                gof = isoterm.gof;
                confInt = confint(isotermResult, confidenceLevel)';
                confInt(confInt<0) = 0;
                this.data = {
                    isotermType.name,...
                    isotermType.formula,...
                    mat2str(coeffvalues(isotermResult), 2),...
                    char(probvalues(isotermResult)),...
                    true,...
                    gof.adjrsquare,...
                    gof.sse,...
                    mat2str(confInt, 2)
                    };
            end
        end
    end
    
end

