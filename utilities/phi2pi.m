function phi = phi2pi(v)
%PHI2PI returns angle phi in [0,2pi] between the vector v=(v_x,v_y) and the x-axis
X=v(1);
Y=v(2);
if X<0 && Y ~=0; phi=atan(Y/X)+pi;
    elseif Y>0 && X>0; phi=atan(Y/X);
    elseif  Y<0 && X>0; phi=atan(Y/X)+2*pi;
    elseif  X==0 && Y>0; phi=pi/2;    
    elseif  X==0 && Y<0; phi=3*pi/2;     
    elseif  X>0 && Y==0; phi=0;
    elseif  X<0 && Y==0; phi=pi;    
end  



end

