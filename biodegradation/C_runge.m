% function y = C_runge(k, T, C0, ssv, v)
%     function res = F(t, y)
%         res = -k(1)*ssv*y/(v*(y+k(2)+y.^2/k(3))); %Haldane
%     end
%
%     options = odeset('RelTol',1e-4, 'AbsTol', 1e-4);
%     [~,y] = ode23t(@F, T, C0, options);
%     y=y';
% end

function y = C_runge(k, T, C0, ssv, v)
    function ress = F(t, y)
        ress = -k(1)*ssv*y/(v*(y+k(2)+y.^2/k(3))); %Haldane
    end
size=length(T);
y=zeros(1,size);
for i=1:size
    if(i == 1)
        y(1)=C0;
    else
        %Real runge kutta
        y(i) = y(i-1) + F(0, y(i-1)) * (T(i)-T(i-1));
        
    end
end
end

