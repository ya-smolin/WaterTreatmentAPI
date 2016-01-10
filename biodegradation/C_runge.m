function y = C_runge(k, T, C0, ssv, v)
    function res = F(t, y)
        res = -k(1)*ssv*y/(v*(y+k(2)+y.^2/k(3))); %Haldane
    end
    
    options = odeset('RelTol',1e-4, 'AbsTol', 1e-4);
    [~,y] = ode23t(@F, T, C0, options);
    y=y';
end

