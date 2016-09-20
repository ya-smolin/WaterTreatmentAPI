function Cout=Cout(t)
global tp CoutM;
i=1;
while(t>=tp(i))&&(i<length(tp))
    i=i+1;
end;
    Cout=CoutM(i);
end
