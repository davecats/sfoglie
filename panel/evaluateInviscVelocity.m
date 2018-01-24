function [u, v] = evaluateInviscVelocity( xi,fld,prf )
%EVALUATEINVISCVELOCITY evaluates the velocity for an arbitrary
%                       Field point


dx=1e-16;
dy=1e-16;

xp=xi;xp(1)=xp(1)+dx;
psix1=evaluateInviscFieldSol(xp,fld,prf);
xm=xi;xm(1)=xm(1)-dx;
psix2=evaluateInviscFieldSol(xm,fld,prf);

u= ( psix2-psix1 )/(2*dx) ;

xp=xi;xp(2)=xp(2)+dy;
psiy1=evaluateInviscFieldSol(xp,fld,prf);
xm=xi;xm(2)=xm(2)-dy;
psiy2=evaluateInviscFieldSol(xm,fld,prf);

v= ( psiy2-psiy1 )/(2*dy) ;

end

