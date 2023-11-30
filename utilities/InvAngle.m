function [ phi ] = InvAngle( theta )
%INVANGLE     Transforms a vectors angle to y-axis into the one to x-axis
%             input:  angle to y-axis theta in [-pi, pi]
%             output: angle to x-axis phi in [-pi, pi] 

scalar=isscalar(theta); % checks wether input is a scalar

if scalar==true
    % third quadrant
    if theta < -pi/2
        phi=-3/2*pi-theta;
    % first, second and fourth quadrant    
    else
        phi=pi/2-theta;
    end
else
    n=length(theta(:,1));
    m=length(theta(1,:));
     % first, second and fourth quadrant
    phi=pi/2*ones(n,m)-theta;
    % third quadrant
    phi(theta < -pi/2)= -3/2*pi -theta(theta < -pi/2); 
end




end