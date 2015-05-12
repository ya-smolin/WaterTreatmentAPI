classdef Model < handle
    
    properties(Constant)
        isotermTypes = {
            %{
            M = mg
            L = litres
            k1: Q0 the maximum sorption density [M/M]
            k2: b the affinity of adsorbent to the adsorbate [L^3/M]
            p1: Sw the aqueous solubility of the sorptive material[M/L]
            %}
            IsotermType('BET','(k1*k2*c)/((p1-c)*(1+(k2-1)*(c/p1)))', [0 1], [Inf, Inf]),...
            %{
            k1: Kf: the Freundlich isotherm parameter[((M/M)/(M/L^3))(1/n)]
            k2: (1/nf): the Freundlich exponent [no units].
            %}
            IsotermType('Freundlich','k1*c^k2', [0 0], [Inf, 1]),...
            %{
            k1: Q0
            k2: b
            k3: Kp: a linear partitioning parameter [L^3/M ]
            %}
            IsotermType('Freundlich-linear','k1*c^k2+k3*c', [0 0 0], [Inf, 1, Inf]),...
            %{
            k1: Q0
            k2: b
            %}
            IsotermType('Langmure','k1*k2*c/(1+k2*c)', [0 0], [Inf, Inf]),...
            %{
            k1: Q0
            k2: b
            k3: (1/ng) : the Generalized Langmuir-Freundlich exponent [no units].
            %}
            IsotermType('Freudlich-Langmure','k1*(k2*c)^k3/(1+(k2*c)^k3)', [0 0 0], [Inf, Inf, 1]),...
            %{
            k1: Q0
            k2: b
            k3: Kp
            %}
            IsotermType('Langmure-Linear','k1*k2*c/(1+k2*c) + k3*c', [0 0 0], [Inf, Inf, Inf]),...
            %{
            k1: Kp
            %}
            IsotermType('Linear','k1*c', [0], [Inf]),...
            %{
            k1: Q0
            k2: A a lumped Polanyi isotherm parameter
            k3: B a lumped Polanyi isotherm parameter
            p1: Sw
            %}
            IsotermType('Polanyi','k1*10^(-k2*(log(p1/c))^k3)', [0 0 0], [Inf, 1, 10]),...
            %{
            k1: Q0
            k2: A
            k3: B
            k4: Kp
            p1: Sw
            %}
            IsotermType('Polanyi-linear','k1*10^(-k2*(log(p1/c))^k3)+k4*c', [0 0 0 0], [Inf, 1, 10, Inf]),...
            %{
            k1: Q0
            k2: b
            k3: nt: the Toth exponent [no units]
            %}
            IsotermType('Toth','k1*k2*c/(1+(k2*c)^k3)^(1/k3)', [0 0 0], [Inf, Inf, 1]),...
            %{
            k1: RT/b
            k2: >=1/min(Cr)
            %}
            IsotermType('Temkin','k1*log(k2*c)', [0 0], [Inf, Inf]),...
            IsotermType('Redlich-Peterson','(k1*k3*c)/(1+k3*c^k2)', [0 0 0], [Inf, 1, Inf]),...
            %{
            k1: theoretical isotherm  saturation  capacity  (mg/g)
            %}
            IsotermType('Dubinina-Radushkevitsa','k1*exp(-k2*log(k3/c)^2)', [0 0 0], [Inf, Inf, Inf]),...
            };
   end
   properties(SetAccess = public)
        Cr;
        Ar;
        parameters;
        lastIsoInd; %crutch for Observable pattern
   end
   
   properties(SetObservable = true)
       isoterms;
   end
   
    methods
        function M = Model(Cr, Ar)
            M.Cr = Cr;
            M.Ar = Ar;
            M.isoterms = cell(length(M.isotermTypes), 1);
        end
        
        function calculate(this, isotermsID)
            for id = isotermsID
                 if id~=1 continue; end
                
                isotermType = this.isotermTypes{id};
                isoterm = Isoterm(isotermType, this.Cr, this.Ar);
                if isempty(isoterm.isotermResult)
                    display(['result is empty for isoterm #' num2str(id)]);
                    continue;
                end
                
                %setter?
                this.lastIsoInd = id;
                this.isoterms{id} = isoterm;
            end
        end
    end
    
end

