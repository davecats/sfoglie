%  
% INPUT PARAMETERS
%

%  profile and panels
%  ------------------
prf.naca = [4 4 1 2];         %  NACA 4-digit profile 
prf.noSkew  = true;           %  if true neglects profile skewness
prf.sharpTE = false;          %  if true modifies NACA prf for sharp trailing edge
prf.c = 1;                    %  prf chord length 
prf.M = 140;                  %  number of nodes on each surface (total number of nodes is 2M-1)
prf.pmode = 1;                %  panelization mode: 1 (more nodes in middle) 2 (more nodes at LE and TE)

%  Flow
%  ----
flo.alfa = 2*pi/180;         %  Angle of attack 
flo.invisc = false;          %  only inviscid solution 
flo.Uinfty=1;                %  velocity of outer flow
flo.Re= 4e5;                 %  Chord Reynoldsnumber
flo.nkrit= 0.15;%            %  critical amplification exponent for transition

%  Tripping
%  --------
tri.active=[ false;...         %  tripping on suction side 
             false];           %  tripping on pressure side 

tri.x = [ 0.145;...            %  tripping location on suction side
          0.29 ]*prf.c;        %  tripping location on pressure side
    
%  Blowing
%  --------
blo.active=false;             % activate blowing
blo.pressureCor=false;        % include correction term for pressure due to blowing -> empirical approximation. Only has small influence.

blo.Mode=1;     % determines how the blowing regions are defined

% blo.Mode=1 -> vector of start x-coordinates and lengths for both sides
% blo.Mode=2 -> vector of mid point x-coordinates and lengths for both sides
% blo.Mode=3 -> Arclength mode 

% Mode 1 and 2
blo.L= {[0.1]*prf.c;              %  length of blowing area: suction side
        [0.1]*prf.c;};            %  length of blowing area: pressure side
blo.x= {[0.4]*prf.c;              %  start/mid point of blowing area
        [0.4]*prf.c;};             
blo.A= {[0.005]*flo.Uinfty;       % blowing intensity of each region 
        [0.005]*flo.Uinfty};

% Mode 3      
% advanced determination of blowing region if blo.Mode==3
%   -> blowing regions are not restricted on upper or lower part
%   -> use intensity 0 to turn off blowing / suction for following region


% commit values of arclength, where blowing intensity changes and the respective velocities
% arclengths in percent of smax -> values have to be lower than 1
%                               -> s_i+1 has to be bigger than s_i
%                               -> no blowing until s_1
blo.s_change=        [0.0233,0.0854,0.458,0.5059]; 
% intensity for s_i< s < s_i+1  -> vector has to be of same length than blo.s_change
blo.NewA= flo.Uinfty*[-0.01,0,0.01,0];
    
    
    
%  Newton Solver
%  --------------
eng.tol=5e-4;                % tolerance of Newton method
eng.it=24;                   % maximum number of Newton step iterations

% Handling of the transition panel   (only for free transition)
% --------------------------------
eng.tranEQ=2;     
% For free transition in few cases some convergence problems occur depending on the mode, that is used
% tranEQ=1 : Xfoil transition EQ 2nd Order. Seems to have convergence problems for very small nkrit and in some blowing cases
% tranEQ=2 : Xfoil transition EQ 1nd Order. Works a bit faster and has less convergence problems but is supposed to be less accurate
% tranEQ=3 : Modified transition EQ. Works often fine but requires more calculation time









