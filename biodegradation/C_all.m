function Y = C_all(k, data)

%     function y = F(k, i)
%         y = -k(1)*ssv*C(i)/(v*(C(i)+k(2)+C(i).^2/k(3))); %Haldane
%         % res = k(1)*X(i)*C(i)/(v*(C(i)+k(2))); %Mono
%     end
 function ress = F(t, y)
        ress = -k(1)*ssv*y/(v*(y+k(2)+y.^2/k(3))); %Haldane
    end
M = length(data);
resG=[];
Call = [];
for testNum = 1:M
    X = data{testNum}.X;  
    T = data{testNum}.T;
    C = data{testNum}.C;
    v = data{testNum}.v;
    Call = [Call C];
    ssv = mean(X);
    
    size=length(T);
    
    res=zeros(1,size);
    for i=1:size
        if(i == 1)
            res(1)=C(1);
        else
            %TODO: review if I correctly take integral of composite fun
            %res(i)= res(i-1)+(F(k, i)+F(k, i-1))./2 .* (T(i)-T(i-1));
             res(i) = res(i-1) + F(0, res(i-1)) * (T(i)-T(i-1));
        end
    end
    
    resG=[resG res];
end

Y = sum((Call - resG).^2);
end

