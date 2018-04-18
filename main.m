init
parameters
[sol,~,~,~,~,~]=airfoil(prf,flo,tri,blo,eng);
CL_ref=sol.CL; CD_ref=sol.Cdrag; tr_ref=sol.tran.x;

N=20;
xTop=linspace(-0.05,1.05,N);
xBot=linspace(-0.05,1.05,N);

tri.x=tr_ref;
tri.active=[true; true];
blo.active=true;

for i=1:N
    for j=1:N
        blo.x= {[xTop(i)]*prf.c;          %  midpoint of blowing area
                [xBot(j)]*prf.c;};  
        [sol,~,~,~,~,~]=airfoil(prf,flo,tri,blo,eng);
        CL(i,j)=sol.CL;
        CD(i,j)=sol.Cdrag;
    end
end

contourf(xTop,xBot,(real(CL(:,:))./CD)'/(CL_ref/CD_ref))
xlabel('x_M^{TOP}')
ylabel('x_M^{BOTTOM}')