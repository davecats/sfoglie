function [ phi ] = InvAngle( theta )
%INVANGLE     Transforms a vectors angle to y-axis into the one to x-axis
%             input:  angle to y-axis theta in [-pi, pi]
%             output: angle to x-axis phi in [-pi, pi] 

scalar=isscalar(theta); % checks wether input is a scalar

if scalar==true
     if theta<-pi/2 %&& theta>-pi
       phi=-3*pi/2-theta;
     elseif theta==0
        phi=0;%pi/2;
     elseif theta==pi/2
        phi=0;
     elseif theta==pi ||theta==-pi
        phi=0;%-pi/2;
     elseif theta==-pi/2
        phi=pi;
    else
       phi=pi/2-theta;
     end
else
    n=length(theta(:,1));
    m=length(theta(1,:));
    phi=pi/2*ones(n,m)-theta;
    phi(theta==-pi)=0;
    phi(theta<-pi/2)=theta(theta<-pi/2)-3*pi/2;
    
    
    phi(theta==0)=0;
    phi(theta==pi/2)=0;
    phi(theta==pi)=0;
    phi(theta==-pi/2)=pi;
    
end













end


