%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%	PROGRAM NAME: ADSORB0.M (MatLab Version 5.2)
% 	SINGLE COMPONENT ADSORPTION KINETICS IN A SINGLE PARTICLE WITH LINEAR ISOTHERM
%  Written by D. D. Do
%             Department of Chemical Engineering
%             University of Queensland, St Lucia, Qld 4072, Australia
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 	           This is the basic problem of solving the problem of adsorption
% 	           kinetics in a single particle with linear isotherm.
%
%	           Ref: D.D.Do,	"Adsorption Analysis: Equilibria and Kinetics"
%			                   Chapter 9. Appendix 9.1
%
%	           Date written: 27 January 1997
%	           Date revised: 15 January 1998:	Convert to MatLab 5.2
%+++++++++++++++++++++++++++++++++++ ++++++++++++++++++++++++++++++++
% NOMENCLATURE:
%
%  Variable            Unit     Description
%  --------            ----     -----------
%	A			         :(-)		: First derivative collocation matrix
%	al			         :(-)		: Parameter alpha for Jacobi polynomial
%	B			         :(-)		: Second derivative collocation matrix
%	be			         :(-)		: Parameter beta for Jacobi polynomial
%	C			         :(-)		: Matrix defined as: 4*u(i)*B(i,j)+2*(s+1)*A(i,j)
%	half_time			:(-)		: Nondimensional half-time for slab, cylinder and sphere
%	n			         :(-)		: Number of interior collocation points
%	n0			         :(-)		: = 1 if the point x=0 is counted; otherwise =0
%	n1			         :(-)		: = 1 if the point x=1 is counted; otherwise =0
%	nt			         :(-)		: n + n0 + n1
%	omega		         :(-)		: Integration vector for ODE15s
%	omega0	         :(-)		: initial integration vector
%	s			         :(-)		: Particle shape factor (0 for slab, 1 for cylinder, 2 for sphere)
%	t, tout	         :(-)		: Nondimensional time
%	t_initial         :(-)		: Initial time (usually = 0)
%	t_final	         :(-)		: Final nondimensional time
%	fractional_uptake	:(-)		: Fractional uptake
%	x			         :(-)		: Nondimensional distance
%	u		            :(-)		: Interpolation points = x^2
%	y			         :(-)		: Nondimensional concentration
%	yb			         :(-)		: Nondimensional bulk concentration
%	w			         :(-)		: Radau quadrature weight
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% USAGE:
%	Type ADSORB0
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% CLEAR ALL VARIABLES
% dbstat=dbstatus;
% save('myBreakpoints.mat', 'dbstat');
% clear all;
% %clear classes;
% close all;
% load('myBreakpoints.mat');
% dbstop(dbstat);
% format shortG;
% clc
    
%pathGen = genpath('C:\WORKSPACE\adsorption book\');
%addpath(pathGen);

% DECLARE SOME GLOBAL VARIABLES
	global A C n u	
	global s
	global yb
	global t_final
   
   %undimensional half time
   %@Ref: page 542, Eq. (9.2-36a)
   half_time   = [0.19674 0.0631 0.03055];%
   
%----------------------------------------------------------------
% USER SUPPLY SECTION
   s           = 2;	                % particle shape factor
   yi          = 0;                 % nondimensional initial concentration
   yb          = 1;                 % nondimensional bulk concentration
   t_initial   = 0;                 % initial time for integration
   t_final     = 20*half_time(s+1); % final time for integration
%----------------------------------------------------------------   
   
% COLLOCATION SECTION
    n			= 10;	
    n0			= 0;
    n1			= 1;
    al          = 1;
    be          = (s-1)/2;

	nt				= n + n0 + n1;
	A				= zeros(nt); 
	B				= zeros(nt); 
	C				= zeros(nt);

	[dif,u] 	   = JCRoot(n,n0,n1,al,be);
    yP = YacobiPolynomial(n, al, be);
    nodes = [yP.u; 1];
    disp ( ['nodes are equal? ' equalEps(u, nodes)]);
    lPs = LagrangePolynomials(nodes);
    derMy=lPs.pNder();
    disp ( ['der are equal? ' equalEps(derMy, dif)]);
    
    Amy=lPs.derValInNodes1(1);
    Bmy=lPs.derValInNodes1(2);
    
	A				= AB(n,n0,n1,1,dif,u);
	B				= AB(n,n0,n1,2,dif,u);
    
    disp ( ['AB are equal? ' equalEps(Amy, A) ' ' equalEps(Bmy, B)]);
    
    w				= RDW(n,n0,n1,al-1,be,u,dif);
    wMy=lPs.rdw(al-1, be-0.5);
    wRight=lPs.rdw00slow();
    testf = @(x)x.^2;
    int1 = wRight'*testf(lPs.nodes);
    %0.0635
    disp ( ['int is as expected? ' equalEps(int1, 1/3)]);
    
    disp ( ['w are equal? ' equalEps(w, wMy)]);
    
	for i=1:nt
		C(i,:)	= 4.*u(i).*B(i,:) + 2.*(s+1).*A(i,:);
	end

% ODE SOLVER
	omega0 			= yi*ones(n,1);
	options			= odeset('Reltol',1e-2,'Abstol',1e-5,'bdf','off');
	[tout,omega] 	= ode15s('Fadsorb0',[t_initial t_final],omega0,options);	

%---------------------------------------------------------------
% SOLUTION: CONCENTRATION PROFILES AND FRACTIONAL UPTAKE
%---------------------------------------------------------------
	yout				= [omega  yb*ones(length(tout),1)];

    
% CALCULATE THE FRACTIONAL UPTAKE
	for i=1:length(tout)
		fractional_uptake(i) 	= ( dot(w, yout(i,:)) - yi )/(yb - yi);
        fractional_uptake1(i) 	= ( dot(wRight, yout(i,:)) - yi )/(yb - yi);
	end

%-----------
% PLOTTING
%-----------
	figure(1)
	subplot(1,2,1)
   plot(tout,fractional_uptake,'k-');
   hold on
   plot(tout,fractional_uptake1,'r-');
   xlabel('time');
   ylabel('Fractional uptake');grid;
   title('FRACTIONAL UPTAKE vs TIME');
   
   subplot(1,2,2)
    dt = floor(length(tout)/10);
    for i=1:dt:length(tout)
		plot(sqrt(u),yout(i,1:n+1),'b-',sqrt(u),yout(i,1:n+1),'bo');
		hold on
   end
   xlabel('Intra-particle distance');
   ylabel('Intra-particle concentration');
   title('CONCENTRATION PROFILES');

   function isEqualStr = equalEps(a, b, epsilon)
       if nargin == 2
           epsilon = eps('single');
       end
       isEqual = sum(sum(abs((a - b) ./ max(a, b)))) < epsilon;
       if(isEqual)
           isEqualStr = 'Yes';
       else
           isEqualStr = 'No';
       end
   end