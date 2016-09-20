function Cin=Cin(t)
global tp CinM;
i=1;
while(t>=tp(i))&&(i<length(tp))
    i=i+1;
end;
    Cin=CinM(i);
end
