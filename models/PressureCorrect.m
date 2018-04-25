function [pInt] = PressureCorrect( s_start,s_end,s,Vb,U ,DT)
%PRESSURECORRECT models additional Pressure Term due to blowing: int dp/dx dy
%                s_start: arclength of blowing start
%                s_end  : arclength of blowing end
%                DT     : Displacement thickness 

dim=size(s);

%ind0=find(s<s_start);
ind1=find(s>=s_start-1e-12 & s<s_end);
ind2=ind1(end)+1:length(s);

DT=DT(ind1);% only in blowing region
pInt=zeros(dim);


if length(Vb)==1
   tmp=Vb;
   Vb=zeros(dim);
   Vb(ind1)=tmp;
end
if length(U)==1
   U=U* ones(dim);
end


% derivate of u*delta_1 for approximation of magnitude of v
% lim_0->infty v-V = d U*delta_1 /ds + v_b
dUDT_ds= FiniteDifferences(U(ind1).*DT,s(ind1));

% normal velocity at BL edge
vBL= dUDT_ds + Vb(ind1);

% linear approximation of v in BL
vmean= 0.5*(vBL + Vb(ind1));

% for first lengthscale
vm1= mean(vmean);
vm2=vmean(end);


Umean=mean(U(ind1));
U(U<0.4*Umean)=0.4*Umean;


%length for effect of Term to expire
k=0.8; % fit to DNS
L1=k*Umean/vm1 *mean(DT); % velocity information reaches displacement thickness
L2=k*U(ind1(end))/vm2 *DT(end);

%sig=-L/log(0.027); -> drop to 3% intensity after L
sig1=L1/3.6;
sig2=L2/3.6;


a1= 2.2*Vb(ind1(1));
a2= 2.2*Vb(ind1(end));

s1=s(ind1)-s_start;
s2=s(ind2)-s_end;

pInt(ind1)= -a1*exp(-s1/sig1) - 0.45*Vb(ind1(1));
pInt(ind2)=  a2*exp(-s2/sig2);
pInt(ind1(end))=0;

%scale
pInt=pInt*0.84;
end

