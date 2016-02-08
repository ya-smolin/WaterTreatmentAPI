function res = C_cal(k, T, C, ssv, v)
    function y = F(k, i)
        global kC
        y = -k(1)*ssv*C(i)/(v*(C(i)+kC(1)+C(i).^2/kC(2))); %Haldane
        % res = k(1)*X(i)*C(i)/(v*(C(i)+k(2))); %Mono
    end
    size=length(T);
    res=zeros(1,size);
    for i=1:size
        if(i == 1)
            res(1)=C(1);
        else
            %TODO: review if I correctly take integral of composite fun
            res(i)= res(i-1)+(F(k, i)+F(k, i-1))./2 .* (T(i)-T(i-1));
        end
    end
end

