function [solnew,Rel, res] = Update(prf,flo,sol,dT, dc,dm,Uinv,D,k)
%UPDATE updates the solution from one Newton step and underrelaxes if necessary

nu=flo.nu;
nkrit=flo.nkrit;


Un=Uinv + D*(sol.m+dm);

% new Lift Coefficient

CLnew = getCL(prf,Un,flo.alfa,flo.Uinfty);
dCL= (CLnew-sol.CL);

dCLmin= max(-0.5, -0.9*sol.CL);

Rel=1;
if dCL > 0.5;  Rel= 0.5/dCL; end
if dCL < dCLmin; Rel=dCLmin/dCL; end

dU= Un  -sol.U;


dDI= (dm-sol.D.*dU)./sol.U; % from product law: m=D*U -> dm=dD*U + D*dU

% rel changes in c
ts=10*ones(size(sol.c)); 
ind=(1: sol.iTran(1)); ts(ind)=sol.c(ind);
ind=(sol.iTran(2):length(sol.c));ts(ind)=sol.c(ind);
DC=dc./ts;
% DC=zeros(size(sol.c)); 
% ind=(1: sol.iTran(1)); DC(ind)=dc(ind)./sol.c(ind);
% ind=(sol.iTran(2):length(sol.c)); DC(ind)=dc(ind)./sol.c(ind);

DT=dT./sol.T;
DD=dDI./sol.D;
DU=dU/0.25;

% residual
% rms= DC.^2 + DT.^2 + DD.^2 + DU.^2;
% res= sqrt( sum(rms)/(4*prf.N+flo.wake.N) );

rms= DC(1:prf.N).^2 + DT(1:prf.N).^2 + DD(1:prf.N).^2 + DU(1:prf.N).^2;
res= sqrt( sum(rms)/(4*prf.N) );

% figure
% hold on
% plot(prf.s,DC(1:prf.N).^2)
% plot(prf.s,DT(1:prf.N).^2)
% plot(prf.s,DD(1:prf.N).^2)
% plot(prf.s,DU(1:prf.N).^2)
% legend('dc','DT','DD','DU')

% adjust Relaxationfactor for big relative changes in D,C,T or U
Rel=GetRelaxationFactor(DC, Rel);
Rel=GetRelaxationFactor(DT, Rel);
Rel=GetRelaxationFactor(DD, Rel);
Rel=GetRelaxationFactor(DU, Rel);


solnew.T=sol.T + Rel*dT;
solnew.D=sol.D + Rel*dDI;
solnew.U=sol.U + Rel*dU;
solnew.c=sol.c + Rel*dc;

solnew.CL=sol.CL + Rel*dCL;

solnew.m=solnew.D.*solnew.U;

%solnew.m=sol.m + Rel*dm; 
%solnew.D=solnew.m./asolnew.U;

% filter unrealistic Ctau values
if sol.iTran(1)~=1    ; turb=(sol.iTran(1)-1:-1:1); else turb=[]; end
turb=[turb,(sol.iTran(2)+1:prf.N+flo.wake.N)];   
tmp=sol.c(turb);
tmp(tmp<0.0000001)=0.0000001;
tmp(tmp>0.25)=0.25;
solnew.c(turb)=tmp;


% filter negativ U values that are not caused by stagnation point movement
%--------------------------------------------------------------------------------------------------
U1= solnew.U(prf.Nle-1:-1:1);
first= find(U1>0,1); % find first positiv velocity
indN=find(U1<0);     % find negative nodes   
indCor=indN(indN>first);% negative nodes that occur after a positiv one have to be corrected
indCor= prf.Nle*ones(size(indCor))-indCor;

U2= solnew.U(prf.Nle:end);
first= find(U2>0,1); % find first positiv velocity
indN=find(U2<0);     % find negative nodes   
indCor2=indN(indN>first);% negative nodes that occur after a positiv one have to be corrected
indCor2= (prf.Nle-1)*ones(size(indCor2))+indCor2;
ind=[indCor; indCor2];

if ~isempty(ind)
    %disp(['UPDATE: negativ U filtered. Nodes ', int2str(ind)])
    sges= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(flo.wake.s)) + flo.wake.s]; 
    for k=1:length(ind)
        if ind(k)==1
            solnew.U(ind(k))=solnew.U(ind(k)+1);
        elseif ind(k)==prf.N+flo.wake.N  || solnew.U(ind(k)+1)<0 %use backward approximation if next velocity is also negativ
            solnew.U(ind(k))=solnew.U(ind(k)-1);
        else % use central approximation if all neighbor velocities are positiv
            solnew.U(ind(k))= ( ( sges(ind(k)+1)-sges(ind(k)) )*solnew.U(ind(k)-1)  +  ( sges(ind(k))-sges(ind(k)-1) )*solnew.U(ind(k)+1)  )...
                                    /( sges(ind(k)+1)-sges(ind(k)-1));
        end
        solnew.m(ind(k))=solnew.U(ind(k))*sol.D(ind(k));
    end
end
%--------------------------------------------------------------------------------------------------



solnew.HK=solnew.D./solnew.T;
solnew.HK(prf.N+1)=(solnew.D(prf.N+1)-prf.gap)./solnew.T(prf.N+1);
solnew.Ret=abs(solnew.U).*solnew.T/nu;
solnew.tran=sol.tran;
solnew.iTran=sol.iTran;
solnew.Vb=sol.Vb;
solnew.Tripping=sol.Tripping;
solnew.sT=sol.sT;
solnew.xT=sol.xT;
% solnew.HmaxLam=sol.HmaxLam;
% solnew.HmaxTurb=sol.HmaxTurb;
solnew.itmax=sol.itmax;
solnew.resmax=sol.resmax;


% eliminate unrealistic values for H
% Hmin=1.02
ind=find(solnew.HK(1:prf.N)<1.02);
dh=max( 0, 1.02*ones(size(ind)) - solnew.D(ind)./solnew.T(ind)  )  ;
solnew.D(ind)=solnew.D(ind) + dh.*solnew.T(ind);
solnew.HK(ind)=solnew.D(ind)./solnew.T(ind);
solnew.m(ind)=solnew.D(ind).*solnew.U(ind);

% on flo.wake: Hmin=1.00005
tmp=solnew.D(prf.N+1:end);tmp(1)=tmp(1) - prf.gap;
HT=tmp./solnew.T(prf.N+1:end);
indw=find(HT<1.00005);
ind=indw + prf.N*ones(size(indw));
dh=max( 0, 1.00005*ones(size(indw)) - tmp(indw)./solnew.T(ind)  );
tmp(indw)=tmp(indw) + dh.*solnew.T(ind);
tmp(1)=tmp(1) +  prf.gap;
solnew.D(ind)=tmp(indw);
solnew.HK(ind)=solnew.D(ind)./solnew.T(ind);
solnew.m(ind)=solnew.D(ind).*solnew.U(ind);

% do not update potential negative values
ind=find(solnew.D<0);
solnew.D(ind)=sol.D(ind);
solnew.HK(ind)=sol.HK(ind);

ind=find(solnew.T<0);
solnew.T(ind)=sol.T(ind);
solnew.HK(ind)=sol.HK(ind);
solnew.Ret(ind)=sol.Ret(ind);

end


















