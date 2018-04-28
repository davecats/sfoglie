%  
% INPUT PARAMETERS
%

%  prf and panels
%  ------------------
prf.naca = [4 4 1 2];         %  NACA 4-digit profile 
prf.noSkew  = true;           %  if true neglects profile skewness
prf.sharpTE = false;          %  if true modifies NACA prf for shart trailing edge
prf.c = 1;                    %  prf chord length 
prf.M = 140;                  %  number of control points for each surface
prf.pmode = 1;                %  panelization mode: 1 (more nodes in middle) 2 (more nodes at LE and TE)

%  Flow
%  ----
flo.alfa = 12*pi/180;         %  Angle of attack 
flo.invisc = false;          %  only inviscid solution 
flo.Uinfty=1;                %  velocity of outer flow
flo.Re= 4e5;                 %  Chord Reynoldsnumber
flo.nkrit= 0.15;%          %  critical amplification exponent for transition

%  Tripping
%  --------
tri.active=[ false;...         %  tripping on suction side 
             false];           %  tripping on pressure side 

tri.x = [ 0.145;...             %  tripping location on suctoin side
          0.29 ]*prf.c;        %  tripping location on pressure side
    
%  Blowing
%  --------
blo.active=false;                 %  activate blowing
blo.L= {[0.1]*prf.c;              %  length of blowing area: suction side
        [0.1]*prf.c;};            %  length of blowing area: pressure side
blo.x= {[0.5]*prf.c;          %  midpoint of blowing area
        [0.5]*prf.c;};             
blo.A= {[0.005]*flo.Uinfty;
        [0.005]*flo.Uinfty};
blo.pressureCor=false;        % include correction term for pressure

%  Newton Solver
%  --------------
eng.tol=5e-4;                % tolerance of Newton method
eng.it=28;                   % maximum number of Newton step iterations
eng.tranEQ=1;                % 1) xfoil transition EQ 2nd Order 2) xfoil transition EQ 1nd Order 3) modified transition EQ
% Transition Equation for free Transition is sometimes a bit tricky when it comes to convergence
% here are some expirience collected with the modes:
%
% tranEQ=1 seems to have convergence problems for small nkrit and in some blowing cases
% tranEQ=2 works a bit faster and has less convergence problems but is supposed to be less accurate
% tranEQ=3 works often fine for small nkrit but not as much for big ones





