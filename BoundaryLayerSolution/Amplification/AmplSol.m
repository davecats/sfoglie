function n2 = AmplSol(flo,n1, T, U, H, Ret, h )
%AMPLSOL solves for the new Amplification of one Interval using a Newton
%method.

% Approximated derivate dn/ds with start value
[dn,~,~,~,~ ]= AmplificationDerivate(flo,H,Ret,T,false );

n2= n1 + max(dn*h,0); % initial guess, make shure there is no decrease


kt=0;res=1;
n=[n1,n2];
% Newton method to solve amplification equation
while res > 1e-9 && kt < 30
    [ f, df_dn,~,~,~, ~ ] = AmplificationEquation(flo,n,T,U,H,Ret,h );

    dn2=-f/df_dn(2);
    % relaxation
    if abs(dn2)>1; rel= 1/abs(dn2); else rel=1; end
    n(2)=n(2) + rel*dn2;
    if n(2)<n1; n(2)=n2; break; end
    
    kt=kt+1;
    res= dn2/n(2);
end
% in case of failed convergence take initial guess
if res<0.1 && n2>n1; n2=n(2); end
      


end

