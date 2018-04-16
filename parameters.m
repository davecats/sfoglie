%  
% INPUT PARAMETERS
%

%  prf and panels
%  ------------------
prf.naca = [4 4 1 2];         %  NACA 4-digit prf 
prf.noSkew = true;           %  if true neglects prf skewness
prf.sharpTE = false ;          %  if true modifies NACA prf for shart trailing edge
prf.c = 1;                    %  prf chord length 
prf.M = 180;                   %  numper of control points for each surface
prf.pmode = 1;                %  panelization mode: 1 (more nodes in middle) 2 (more nodes at LE and TE)

%  Flow
%  ----
flo.alfa = 0*pi/180;         %  Angle of attack 
flo.invisc = false;          %  only inviscid solution 
flo.Uinfty=1;                %  velocity of outer flow
flo.Re= 4e5;                 %  Chord Reynoldsnumber
flo.nkrit=  9;%           %  critical amplification exponent for transition

%  Tripping
%  --------
tri.active=[ false;...         %  tripping on suction side 
             false];           %  tripping on pressure side 

tri.x = [ 0.07;...             %  tripping location on suctoin side
          0.19 ]*prf.c;        %  tripping location on pressure side
    
%  Blowing
%  --------
blo.active=[false;...    %  blowing on suction side  
            false];      %  blowing on pressure side 

blo.xStart= [0.25;...         %  start of blowing on suction side
             0.25]* prf.c;    %  start of blowing on pressure side
   
blo.xEnd  = [0.5;...          %  end of blowing on suction side
             0.5]* prf.c;     %  end of blowing on pressure side
          
blo.A=[0.001;...              %  blowing intensity suction side
       0.001]* flo.Uinfty;    %  blowing intensity pressure side

blo.pressureCor=false;        % include correction term for pressure

%  Newton Solver
%  --------------
eng.it=25;                   % maximum number of Newton step iterations
eng.tranEQ=false;            % false) xfoil transition EW, true) modified transition EQ
eng.tol=5e-4;                % tolerance of Newton method
 
