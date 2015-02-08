classdef Model < handle
    
    properties
        isotermTypes = {IsotermType('Langmure','k1*k2*c/(1+k2*c)', [0 0], [Inf, Inf]),...
            IsotermType('Freundlich','k1*c^k2', [0 0], [Inf, 1]),...
            IsotermType('Freudlich-Langmure','k1*(k3*c)^k2/(1+(k3*c)^k2)', [0 0 0], [Inf, 1, Inf]),...
            IsotermType('Temkin','k1*log(k2*c)', [0 0], [Inf, Inf]),...
            IsotermType('Redlich-Peterson','(k1*k3*c)/(1+k3*c^k2)', [0 0 0], [Inf, 1, Inf]),...
            IsotermType('Dubinina-Radushkevitsa','k1*exp(-k2*log(k3/c)^2)', [0 0 0], [Inf, Inf, Inf]),...
            IsotermType('BET','(k1*k2*c)/((1-c/p1)*(1+(k2-1)*(c/p1)))', [0 0], [Inf, Inf])};
        Cr;
        Ar;
        isoterms;
   end
    
    methods
        function M = Model(Cr, Ar)
            M.Cr = Cr;
            M.Ar = Ar;
            M.isoterms = cell(length(M.isotermTypes), 1);
        end
        
        function calculate(M)
            %global output;
            %output = [];
            %sum = 0;
            M.isotermTypes{7}.constants = 56;
            for i = 1:length(M.isotermTypes)
                if i~=7 continue; end
                %sum = sum + gof.sse;
                M.isoterms{i} = IsotermModel(M.isotermTypes{i}, M.Cr, M.Ar);
                if isempty(M.isoterms{i}.isotermResult)
                    continue;
                end
                subplot(3,3,i)
                plot(M.Cr, M.Ar,'greeno', 'LineWidth', 3);
                hold on;
                plot(M.isoterms{i}.isotermResult)      
            end
            %disp(sum);   
        end
    end
    
end

