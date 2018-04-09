function [pInt] = PressureCorrect( s_start,s_end,s,Vb,U )
%PRESSURECORRECT 

dim=size(s);



%ind0=find(s<s_start);
ind1=find(s>s_start & s<s_end);
ind2=find(s>s_end);

pInt=zeros(dim);


if length(Vb)==1
   tmp=Vb;
   Vb=zeros(dim);
   Vb(ind1)=tmp;
end
if length(U)==1
   U=U* ones(dim);
end

Vb(ind2)=max(Vb);

delta=s_end-s_start;
%sig=-delta/log(0.01);
sig=delta/3.6;

a= 2.2*Vb./U;

s1=s(ind1)-s_start;
s2=s(ind2)-s_end;

pInt(ind1)= -a(ind1).*exp(-s1/sig);
pInt(ind2)=  a(ind2).*exp(-s2/sig);

end

