function res = C_cal(k, T)
    global C
    size=length(T);
    res=zeros(1,size);
    for i=1:size
        if(i == 1)
            res(1)=C(1);
        else
            res(i)= res(i-1)-(F(k, i)+F(k, i-1))./2 .* (T(i)-T(i-1));
        end
    end
end

function res = F(k, i)
global X C v
%res = k(3)*k(1)*X(i)*C(i)/(v*(k(3)*C(i)+k(3)*k(2)+C(i).^2));
res = k(1)*X(i)*C(i)/(v*(C(i)+k(2)+C(i).^2/k(3))); %Haldane
% res = k(1)*X(i)*C(i)/(v*(C(i)+k(2))); %Mono
end