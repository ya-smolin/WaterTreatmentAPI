function W = polinterp(xcol,xint)
  ncol = length(xcol);
  nint = length(xint);
  W = zeros(nint, ncol);
  for i = 1: ncol
    d = 1.0;
    for j = 1: ncol
      if (j ~= i )
	d = d*(xcol(i) - xcol(j));
      end
    end
    for k = 1: nint
      xx = xint(k);
      n = 1.0;
      for j = 1: ncol
	if (j ~= i)
	  n = n*(xx- xcol(j));
	end
      end
      W(k,i) = n/d;
    end
  end  
