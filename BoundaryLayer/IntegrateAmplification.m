function [ n2 ] = IntegrateAmplification(n1,h,T1,T2,H1,H2,Ret1,Ret2)
%INTEGRATEAMPLIFICATION integrates the amplification factor for node 2



H=[H1;H2];T=[T1;T2]; Ret=[Ret1;Ret2];

[dn ,~,~,~] = AmplificationDerivate(H,Ret,T );


n2=n1 + dn*h;


end

