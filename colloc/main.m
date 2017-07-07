%% Copyright (C) 2001, James B. Rawlings and John G. Ekerdt
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation; either version 2, or (at
%% your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, write to
%% the Free Software Foundation, 59 Temple Place - Suite 330, Boston,
%% MA 02111-1307, USA.

global Phi A B R n

Phi  = 10;
n    = 1;
nint = 100;
xint = linspace(0, 3, nint)';


colvec = [5; 10; 30; 50];
ncol   = length(colvec);
ccol   = zeros(ncol, ncol);
csave  = zeros(nint, ncol);
cerr   = zeros(nint, ncol);
can    = 3./xint.*sinh(Phi*xint)/sinh(3*Phi);
can(1) = 3*Phi/sinh(3*Phi);
etaan = 1./Phi*(1./tanh(3*Phi)-1/(3*Phi));
for i = 1: ncol
  npts = colvec(i);
  %%
  %% global collocation 
  %%
  [R A B Q] = colloc(npts-2, 'left', 'right');
  R = R*3;
  A = A/3;
  B = B/9;
  Q = Q*3;
  Aint = A(2:npts-1,:);
  Bint = B(2:npts-1,:);
  Rint = R(2:npts-1);
  W = polinterp(R,xint);
  %%
  %% solve the problem
  %%
  c0=0.5*linspace(0,1,npts)';
  tol = 1e-10;
  opts = optimset ('TolFun', tol);
  [c,fval,info] = fsolve('pellet',c0,opts);
  cint = W*c;
  csave(:,i) = cint;
  cerr(:,i)  = abs((cint - can))./can;
  %%
  %% compute the effectiveness factor
  %% 
  eta  = (n+1)/2*A(npts,:)*c/Phi^2;
  etaerr(i) = abs(eta - etaan)/etaan;
  results = [R c];
%  save ('-ascii', sprintf ('pelletcolloc_%d.dat', i), 'results');
end
table = [xint csave cerr];
save -ascii pelletcolloc.dat table;

%%plot the pellet profiles for different collocation numbers
plot(xint, csave(:,1:ncol))
title ('Figure A.5')
