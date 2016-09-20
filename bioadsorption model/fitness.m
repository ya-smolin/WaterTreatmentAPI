function fit = fitness(fun,x,y)
N=length(x);
fitness=zeros(1,N);
for i=1:N;
fitness(i)=abs(fun(x(i))-y(i));
fit=sum(fitness);
end

