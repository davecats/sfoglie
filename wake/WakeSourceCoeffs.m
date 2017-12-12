function [dpsi_dq,du_dq] = WakeSourceCoeffs( wake , xi, n)
%WAKESOURCECOEFFS  calculates coefficients of contribution of the wake source
%                    -on the airfoil B_ij=dpsi_dq
%                    -on the airfoil Cq_ij=du_dq
%                  source is linear for the wake !!!!

NW=length(wake.x);% number of nodes
L=wake.L;
t= wake.e; %tangential vectors

dpsi_dq=zeros(1,NW);
du_dq=zeros(1,NW);

for k=1:NW-1

    r1=[xi(1);xi(2)]-[wake.x(k);wake.y(k)];% relativ vector

    % calculate local coordinates
    Y=wake.n(k,:)*r1;
    X1=wake.e(k,:)*r1;
    XM=X1-L(k)/2;
    X2=X1-L(k);

    rM=sqrt(XM^2+Y^2);
    r1=sqrt(X1^2+Y^2);
    r2=sqrt(X2^2+Y^2);

     sgnM=sign(XM);sgn1=sign(X1); sgn2=sign(X2);
     if rM<1e-10; lnrM=0; tM=0;  else lnrM=log(rM);  tM=atan2(sgnM*XM,sgnM*Y)+pi/2*(1-sgnM); end
     if r1<1e-10; lnr1=0; t1=0;  else lnr1=log(r1);  t1=atan2(sgn1*X1,sgn1*Y)+pi/2*(1-sgn1); end
     if r2<1e-10; lnr2=0; t2=0;  else lnr2=log(r2);  t2=atan2(sgn2*X2,sgn2*Y)+pi/2*(1-sgn2); end
    tM=tM- InvAngle(atan2(-t(1),t(2) )+pi);
    t1=t1- InvAngle(atan2(-t(1),t(2) )+pi);
    t2=t2- InvAngle(atan2(-t(1),t(2) )+pi);

    % first panel half
    %--------------------------------------------------
    
    % dpsi_dq
    pp=(XM*tM-X1*t1+Y*(lnr1-lnrM));
    pm=((X1+XM)*pp+r1^2*t1-rM^2*tM+Y*(XM-X1))/(X1-XM);

    if k==1
    a=1/L(1);b=a;
    else
    a=1/(L(k)+L(k-1));
    b=1/L(k);
    end
    c1=(-pp*a+pm*a)/(4*pi);
    c2=(-pp*b-pm*b)/(4*pi);
    c3=(pp*(b+a)+pm*(b-a))/(4*pi);
    
    if k==1
        dpsi_dq(k:k+1)=dpsi_dq(k:k+1)+[c1+c2, c3];
    else
        dpsi_dq(k-1:k+1)=dpsi_dq(k-1:k+1)+[c1, c2, c3];
    end
    
    % du_dq
    if nargin==3
        xN= sum(wake.e(k,:).*n);
        yN= wake.e(k,1)*n(2)-wake.e(k,2)*n(1);

        up=-t1*xN+tM*xN + (lnr1-lnrM)*yN;
        tmp1=((X1-XM)*t1 + pp  - pm)/(X1-XM);
        tmp2=((X1-XM)*tM + pp  + pm)/(X1-XM);
        tmp3=((X1+XM)*(lnr1-lnrM)+2*(X1-XM+Y*(t1-tM)));
        um=-tmp1*xN+tmp2*xN + (lnr1-lnrM)*tmp3;
        c1=(-up*a+um*a)/(4*pi);
        c2=(-up*b-um*b)/(4*pi);
        c3=(up*(b+a)+um*(b-a))/(4*pi);
        if k==1
            du_dq(k:k+1)=dpsi_dq(k:k+1)+[c1+c2, c3];
        else
            du_dq(k-1:k+1)=dpsi_dq(k-1:k+1)+[c1, c2, c3];
        end
    end
    
    % second panel half
    %--------------------------------------------------
    pp=(X2*t2-XM*tM+Y*(lnrM-lnr2));
    pm=((XM+X2)*pp+rM^2*tM-r2^2*t2+Y*(X2-XM))/(XM-X2);

    if k==NW-1
    c=b;
    else
    c=1/(L(k+1)+L(k));
    end
    c1=(-pp*(c+b)-pm*(c-b))/(4*pi);
    c2=(pp*b-pm*b)/(4*pi);
    c3=(pp*c+pm*c)/(4*pi);
    if k==NW-1
        dpsi_dq(k:k+1)=dpsi_dq(k:k+1)+[c1, c2+c3];
    else
        dpsi_dq(k:k+2)=dpsi_dq(k:k+2)+[c1, c2, c3];
    end

    if nargin==3
        up=-tM*xN+tM*xN + (lnrM-lnr2)*yN;
        tmp1=((X2-XM)*tM + pp  - pm)/(XM-X2);
        tmp2=((X2-XM)*t2 + pp  + pm)/(XM-X2);
        tmp3=((X2+XM)*(lnrM-lnr2)+2*(X2-XM+Y*(tM-t2)));
        um=-tmp1*xN+tmp2*xN + (lnrM-lnr2)*tmp3;
        c1=(-up*a+um*a)/(4*pi);
        c2=(-up*b-um*b)/(4*pi);
        c3=(up*(b+a)+um*(b-a))/(4*pi);
        if k==NW-1
            du_dq(k:k+1)=du_dq(k:k+1)+[c1, c2+c3];
        else
            du_dq(k:k+2)=du_dq(k:k+2)+[c1, c2, c3];
        end
    end
    
end

end

