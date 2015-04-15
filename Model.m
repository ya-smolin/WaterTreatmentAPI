classdef Model < handle
    
    properties
        isotermTypes = {
            %{
            k1: Q0 the maximum sorption density [M/M]
            k2: b the affinity of adsorbent to the adsorbate [L^3/M]
            p1: Sw the aqueous solubility of the sorptive material
            %}
            IsotermType('BET','(k1*k2*c)/((p1-c)*(1+(k2-1)*(c/p1)))', [0 1], [1e6, 1e6]),...
            %{
            k1: Kf: the Freundlich isotherm parameter[((M/M)/(M/L^3))(1/n)]
            k2: (1/nf): the Freundlich exponent [no units].
            %}
            IsotermType('Freundlich','k1*c^k2', [0 0], [1e6, 1]),...
            %{
            k1: Q0
            k2: b
            k3: Kp: a linear partitioning parameter [L^3/M ]
            %}
            IsotermType('Freundlich-linear','k1*c^k2+k3*c', [0 0 0], [1e6, 1, 1e6]),...
            %{
            k1: Q0
            k2: b
            %}
            IsotermType('Langmure','k1*k2*c/(1+k2*c)', [0 0], [1e6, 1e6]),...
            %{
            k1: Q0
            k2: b
            k3: (1/ng) : the Generalized Langmuir-Freundlich exponent [no units].
            %}
            IsotermType('Freudlich-Langmure','k1*(k2*c)^k3/(1+(k2*c)^k3)', [0 0 0], [1e6, 1e6, 1]),...
            %{
            k1: Q0
            k2: b
            k3: Kp
            %}
            IsotermType('Langmure-Linear','k1*k2*c/(1+k2*c) + k3*c', [0 0 0], [1e6, 1e6, 1e6]),...
            %{
            k1: Kp
            %}
            IsotermType('Linear','k1*c', [0], [1e6]),...
            %{
            k1: Q0
            k2: A a lumped Polanyi isotherm parameter
            k3: B a lumped Polanyi isotherm parameter
            p1: Sw
            %}
            IsotermType('Polanyi','k1*10^(-k2*(log(p1/c))^k3)', [0 0 0], [1e6, 1, 10]),...
            %{
            k1: Q0
            k2: A
            k3: B
            k4: Kp
            p1: Sw
            %}
            IsotermType('Polanyi-linear','k1*10^(-k2*(log(p1/c))^k3)+k4*c', [0 0 0 0], [1e6, 1, 10, 1e6]),...
            %{
            k1: Q0
            k2: b
            k3: nt: the Toth exponent [no units]
            %}
            IsotermType('Toth','k1*k2*c/(1+(k2*c)^k3)^(1/k3)', [0 0 0], [1e6, 1e6, 1]),...
            %{
            k1: RT/b
            k2: >=1/min(Cr)
            %}
            IsotermType('Temkin','k1*log(k2*c)', [0 0], [1e6, 1e6]),...
            IsotermType('Redlich-Peterson','(k1*k3*c)/(1+k3*c^k2)', [0 0 0], [1e6, 1, 1e6]),...
            %{
            k1: theoretical isotherm  saturation  capacity  (mg/g)
            %}
            IsotermType('Dubinina-Radushkevitsa','k1*exp(-k2*log(k3/c)^2)', [0 0 0], [1e6, 1e6, 1e6]),...
            };
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
            %M.isotermTypes{7}.constants = 56;
            
            for i = 1:length(M.isotermTypes)
%                 if i~=1 continue; end  
                isotermType = M.isotermTypes{i};
                curIsotermFitmodel = isotermType.fitmodel;
                curConstants = isotermType.constants;
                pNum = length(probnames(curIsotermFitmodel));
                if(isempty(curConstants) && pNum~=0)
                    disp('avto parameter fitting');
                    fstr = formula(isotermType.fitmodel);
                    kNum = numcoeffs(curIsotermFitmodel);
                    lowerB = isotermType.lowerB;
                    upperB = isotermType.upperB;
                    for j = 1:pNum
                        fstr = strrep(fstr, strcat('p',num2str(j)), strcat('k',num2str(kNum+1)));
                        lowerB = [lowerB, max(M.Cr)];
                        upperB = [upperB, Inf];
                    end
                    isotermType = IsotermType(isotermType.name, fstr, lowerB, upperB);
                end
                M.isotermTypes{i} = isotermType; %change modelType or not?
                if(strcmp(isotermType.name,'Temkin'))
                    isotermType.lowerB = [0, 1/min(M.Cr(2:end))];
                    M.isoterms{i} = IsotermModel(isotermType, M.Cr(2:end), M.Ar(2:end));
                else
                     M.isoterms{i} = IsotermModel(isotermType, M.Cr, M.Ar);
                end
               
                if isempty(M.isoterms{i}.isotermResult)
                    continue;
                end
                
                subplot(3,3,mod(i-1, 9)+1)
                plot(M.Cr, M.Ar,'greeno', 'LineWidth', 3);
                hold on;
                plot(M.isoterms{i}.isotermResult);
                display(M.isoterms{i}.isotermResult);
                legend('off');
            end
        end
    end
    
end

