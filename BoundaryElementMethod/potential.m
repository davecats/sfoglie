function [flo] = potential(prf,flo)
%POTENTIAL  calculates the coefficient matrix A of the equation system 
%           A_ij z_j = -U y_i + V x_i with solution vector z=[gamma_1,..gamma_M, psi_0]
%           and solves it

% right hand side
t=transpose(-flo.ui*prf.nodes.Y + flo.vi*prf.nodes.X);
flo.t=t;
N=length(t)-1;


% coefficient matrix
[A,B]=VortexAndSourceDistr(prf);
flo.A=A; flo.B=B;

%Add streamfunction psi0 as unknown
A=[A, -ones(N+1,1)];

if prf.sharpTE   
     A(end,:)=zeros(1,N+2); % delete Equation from last node -> equals to first node
     t(end)=0;
     
     % no reverse flow condition
     A(end,2)=-2;A(end,3)=1;A(end,end-2)=2;A(end,end-3)=-1; 
     A(end,1)=1;A(end,end-1)=-1;
end
% add Kutta Condition
t=[t;0];
A=[A; zeros(1,N+2)];
A(end,[1, N+1])=[1,1];


Lsg=A\t;    %solve the equation system
flo.gamma = Lsg(1:end-1); % gammas -> solution vector without psi0 -> last element
flo.psi0=Lsg(end);
flo.Ages=A;

end

