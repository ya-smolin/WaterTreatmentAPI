function retval = pellet(c)
  global Phi A B R n
  retval = B*c + 2*A*c./R - 2/(n+1)*Phi^2*(c.^n);
  retval(1)    = A(1,:)*c;
  retval(end) = 1 - c(end);
