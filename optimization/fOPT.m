function  ratio  = fOPT( xB)
%FOPT function for optimization

global OptiInputStruct

blo=OptiInputStruct.blo;
% overwrite parameters with blowing configuration

blo.active=true;                 %  activate blowing
blo.L= {[0.1]*OptiInputStruct.prf.c;              %  length of blowing area: suction side
        [0.1]*OptiInputStruct.prf.c;};            %  length of blowing area: pressure side
blo.x= {xB(1)*OptiInputStruct.prf.c;          %  midpoint of blowing area
        xB(2)*OptiInputStruct.prf.c;};             
blo.A= {[0.005]*OptiInputStruct.flo.Uinfty;
        [0.005]*OptiInputStruct.flo.Uinfty};
blo.pressureCor=false;        % include correction term for pressure


[sol,~,~,~,~,~]=airfoil(OptiInputStruct.prf,OptiInputStruct.flo,...
                        OptiInputStruct.tri,blo,OptiInputStruct.eng,OptiInputStruct);

%disp(['Upper ',num2str(xB(1)),' lower ',num2str(xB(2)), ' ratio ',num2str(sol.CL/sol.Cdrag)])

ratio=-sol.CL/sol.Cdrag;

%make sure unconvereged solution is not found as solution
if sol.residual>1e-3;
   ratio=0; 
end

end

%PlotProfile(prf,flo.wake,sol, 3);

