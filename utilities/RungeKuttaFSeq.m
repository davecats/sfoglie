function [ lsg , diverged] = RungeKuttaFSeq( beta, Initial, h, Inv )
%RUNGEKUTTAFSeq Integriert die Falkner-Skan-Gleichung mit einem Runge-Kutta
%               Verfahren (für Schießverfahren benötigt)
%               Input:  -beta, 
%                       -Anfangswerte für f, f', f''
%                       -Vektor mit Schrittweiten
%                       -Inv falls umgekehrte Keilströmung


if nargin==3;
    Inv=false;
end

sI=size(Initial);
if sI(2)>1
    Initial=Initial';
end
sh=size(h);
if sh(1)>1
    h=h';
end
N=length(h)+1;

lsg= zeros(length(Initial), N);
lsg(:,1)=Initial;

phi1=Initial;

k1=zeros(size(phi1));
k2=k1;k3=k1;k4=k1;

diverged=false;

if Inv % falls umgekehrte Keilströmung -> Vorzeichen anpassen
    sgn=-1;
else
    sgn=1;
end


for i=2:N
    % 1. Faktor -> Gleichung an Stelle i
    k1(1)= phi1(2);
    k1(2)= phi1(3);
    k1(3)=-sgn*phi1(1)*phi1(3) - sgn*beta*(1-phi1(2)^2);
    
    % Eulersches polygonzugverfahren für h/2
    phiTMP= phi1 + h(i-1)/2 * k1;
    
    k2(1) = phiTMP(2);
    k2(2) = phiTMP(3);
    k2(3) = -sgn*phiTMP(1)*phiTMP(3) - sgn*beta*(1-phiTMP(2)^2);
    
    phiTMP= phi1 + h(i-1)/2 * k2;
    
    k3(1) = phiTMP(2);
    k3(2) = phiTMP(3);
    k3(3) = -sgn*phiTMP(1)*phiTMP(3) - sgn*beta*(1-phiTMP(2)^2);
    
    phiTMP= phi1 + h(i-1) * k3;
    
    k4(1) = phiTMP(2);
    k4(2) = phiTMP(3);
    k4(3) = -sgn*phiTMP(1)*phiTMP(3) - sgn*beta*(1-phiTMP(2)^2);
    
    
    phi2= phi1 + h(i-1)/6 * (k1 + 2*k2 + 2*k3 + k4);
    lsg(:,i)=phi2;
    
    % neues phi1 für nächsten Iterationsschritt
    phi1=phi2;
    
    % bricht ab falls die Lösung divergiert
    if abs(phi1(2))>100 
       diverged=true; 
       break; 
    end
    
end





end

