function [n, TranU,TranL] = IntegrateAmplification(prf,sol,flo)
%INTEGRATEAMPLIFICATION     integrates the amplification factor for each
%boundary node with known shape parameter and Ret values

nkrit=flo.nkrit;

% suction side

n=zeros(size(sol.c));
ind=(prf.Nle-1:-1:sol.TranU);
indM=ind(2:end);

dn = AmplificationDerivate(flo,sol.HK(ind),sol.Ret(ind),sol.T(ind) );
dn=max(0,dn); % make sure that Amplification does not decrease

h=prf.panels.L(indM)';

n(indM)=cumsum(dn.*h);

TranU=max(find(n>nkrit,1,'last'));

%pressure side
ind=(prf.Nle:sol.TranL);
indM=ind(2:end);

dn = AmplificationDerivate(sol.HK(ind),sol.Ret(ind),sol.T(ind) );

h=prf.panels.L(ind(1:end-1))';

n(indM)=cumsum(dn.*h);

np=n(prf.Nle:prf.N);
TranL=find(np > nkrit,1);
TranL=TranL+prf.Nle-1;

if isempty(TranU); TranU=1;end
if isempty(TranL); TranL=prf.N;end

end

