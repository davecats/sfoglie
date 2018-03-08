function [ Ts, Ds ] = ImproveStartNodeSol( U,L,Vb, T, D )
%IMPROVESTARTNODESOL improves the start node solution


if nargin==3 % Eppler start values
    nu=evalin('base','nu');
    T= 0.29004*sqrt(L*nu/U); 
    D= 2.0754*T;
end

Ts=T;
Ds=D;

k=0;res=1;
while res>1e-9 && k<60
    [f1,f2, J ] = InitialNodeSys(Ts,Ds,U,Vb,L);
    J=J(:,1:2);
    
    % correction of old solution x(k+1)=x(k) + dx(k)
    % solution of the Equation J dx(k)=-f(k) with cramers law  
    dT= (f2*J(1,2) - f1*J(2,2))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );
    dD= (f1*J(2,1) - f2*J(1,1))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );
    
    res=max(abs( [dT/Ts,dD/Ds] ));
    % under relaxation for big changes
    if res>0.3; Rel=0.3/res; else Rel=1; end            
    % update values
    Ts=Ts + Rel*dT;
    Ds=Ds + Rel*dD;
    
    k=k+1;
end

if res>0.05 %  no convergence take old values
    Ts=T;
    Ds=D;
end



end

