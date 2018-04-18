%  
% INPUT PARAMETERS
%

%  prf and panels
%  ------------------
prf.naca = [4 4 1 2];         %  NACA 4-digit prf 
prf.noSkew  = false;           %  if true neglects prf skewness
prf.sharpTE = false;          %  if true modifies NACA prf for shart trailing edge
prf.c = 1;                    %  prf chord length 
prf.M = 90;                  %  number of control points for each surface
prf.pmode = 1;                %  panelization mode: 1 (more nodes in middle) 2 (more nodes at LE and TE)

%  Flow
%  ----
flo.alfa = 5*pi/180;         %  Angle of attack 
flo.invisc = false;          %  only inviscid solution 
flo.Uinfty=1;                %  velocity of outer flow
flo.Re= 4e5;                 %  Chord Reynoldsnumber
flo.nkrit=  9;%          %  critical amplification exponent for transition

%  Tripping
%  --------
tri.active=[ true;...         %  tripping on suction side 
             true];           %  tripping on pressure side 

tri.x = [ 0.145;...             %  tripping location on suctoin side
          0.29 ]*prf.c;        %  tripping location on pressure side
    
%  Blowing
%  --------
blo.active=false;                 %  activate blowing
blo.L= {[0.1]*prf.c;              %  length of blowing area: suction side
        [0.1]*prf.c;};            %  length of blowing area: pressure side
blo.x= {[0.4687]*prf.c;          %  midpoint of blowing area
        [0.1233]*prf.c;};             
blo.A= {[0.005]*flo.Uinfty;
        [0.005]*flo.Uinfty};
blo.pressureCor=false;        % include correction term for pressure

%  Newton Solver
%  --------------
eng.it=20;                   % maximum number of Newton step iterations
eng.tranEQ=false;            % false) xfoil transition EW, true) modified transition EQ
eng.tol=5e-4;                % tolerance of Newton method
 
