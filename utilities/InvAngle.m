function [ phi ] = InvAngle( theta)
%INVANGLE     Transforms a vectors angle to y-axis into the one to x-axis
%             input:  angle to y-axis theta in [-pi, pi]
%             output: angle to x-axis phi in [-pi, pi] 



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


% if y<0 && x<0
%    phi=-3*pi/2-atan2(y,x);
% elseif x>0 && y==0
%    phi=pi/2;
% elseif x==0 && y>0
%    phi=0;
% elseif x<0 && y==0
%    phi=pi;
% elseif x==0 && y<0
%    phi=-pi/2;
% else
%    phi=pi/2-atan2(y,x); 
% end




end


