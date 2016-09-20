function v=v(t)
global tp vp;
i=1;
while((t>=tp(i))&&(i<length(tp)))
    i=i+1;
end;
    v=vp(i);
end
