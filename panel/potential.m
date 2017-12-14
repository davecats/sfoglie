function [ fld ] = potential( prf)
%POTENTIAL  calculates the coefficient matrix A of the equation system 
%           A_ij z_j = -U y_i + V x_i with solution vector z=[gamma_1,..gamma_M, psi_0]
%           and solves it

% right hand side
u=evalin('base','ui');
v= evalin('base','vi');

t=transpose(-u*prf.nodes.Y+v*prf.nodes.X);
fld.t=t;
N=length(t)-1;


% coefficient matrix
[A,B]=VortexAndSourceDistr(prf);

fld.A=A; 
fld.B=[B(:,1:1:end-1), zeros(N+1,1)]; % linear ansatz, last node doesn´t contribute 
%                                     -> q_TE depends on gamma_1 and is part of A


%Add streamfunction psi0 as unknown
A=[A, -ones(N+1,1)];

% add Kutta Condition
t=[t;0];
A=[A; zeros(1,N+2)];
A(end,[1, N+1])=[1,1];


Lsg=A\t;    %solve the equation system
fld.gamma = Lsg(1:N+1); % gammas -> solution vector without psi0 -> last element
fld.psi0=Lsg(N+2);


end

