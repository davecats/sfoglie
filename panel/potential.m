function [ fld ] = potential( prf)
%POTENTIAL  calculates the coefficient matrix A of the equation system 
%           A_ij z_j = -U y_i + V x_i with solution vector z=[gamma_1,..gamma_M, psi_0]
%           and solves it

% right hand side
u=evalin('base','ui');
v= evalin('base','vi');

t=transpose(-u*prf.nodes.Y + v*prf.nodes.X);
fld.t=t;
N=length(t)-1;


% coefficient matrix
[A,B]=VortexAndSourceDistr(prf);

fld.A=A; 
fld.B=B;

%Add streamfunction psi0 as unknown
A=[A, -ones(N+1,1)];

if prf.IsSharp
%     arg1=atan2(-prf.nodes.Y(1),-prf.nodes.X(1));
%     arg2=atan2(prf.nodes.Y(1),prf.nodes.X(1));
%     arg2=arg2-2*pi*round(arg2-arg1+0.5*sign(arg2-arg1));
%     arg=(arg1+arg2)/2;
%     Cb=cos(arg);Sb=sin(arg);
%     xB=prf.nodes.X(1) - 0.1*min(prf.panels.L(1),prf.panels.L(2))*Cb;
%     yB=prf.nodes.Y(1) - 0.1*min(prf.panels.L(1),prf.panels.L(2))*Sb; 
    
     A(end,:)=zeros(1,N+2); % delete Equation from last node -> equals to first node
     t(end)=0;
     A(end,2)=-2;A(end,3)=1;A(end,end-2)=2;A(end,end-3)=-1; 
     A(end,1)=1;A(end,end-1)=-1;
end
% add Kutta Condition
t=[t;0];
A=[A; zeros(1,N+2)];
A(end,[1, N+1])=[1,1];


Lsg=A\t;    %solve the equation system
fld.gamma = Lsg(1:end-1); % gammas -> solution vector without psi0 -> last element
fld.psi0=Lsg(end);
%if prf.IsSharp; fld.gamma=[fld.gamma ;Lsg(1)];

fld.Ages=A;

end

