function phi = phipi(v)
%PHIPI returns angle phi in [-pi,pi] between the vector v=(v_x,v_y) and the x-axis
X=v(1);
Y=v(2);
if Y>0 && X>0; phi=atan(Y/X);
    elseif X<0 && Y>0; phi=atan(Y/X)+pi;
    elseif X<0 && Y<0; phi=atan(Y/X)-pi;    
    elseif  X>0 && Y<0; phi=atan(Y/X);
    elseif  X==0 && Y>0; phi=pi/2;    
    elseif  X==0 && Y<0; phi=-pi/2;     
    elseif  X>0 && Y==0; phi=0;
    elseif  X<0 && Y==0; phi=pi;    
end  



end

