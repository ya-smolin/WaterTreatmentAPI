function y = C_all(k, T, C0, ssv, v)
    function ress = F(t, y)
         global kC
        ress = -k(1)*ssv*y/(v*(y+kC(1)+y.^2/kC(2))); %Haldane
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
