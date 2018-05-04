function  ratio  = fOPT2(z)
%FOPT2 function for optimization of suction side blowing.
%      Takes start point, length and Intensity as optimization parameters


% the algorithm can not deal with variable constraints so that kind of dirty
% trick has to be used. setting ratio to 0 will make the algorithm go away
% from solutions where blowing end point would be out of the bounds
if z(1)+z(2)>1 
    ratio=0; 
    return;
end


% calculate midpoint
xm=z(1)+z(2)/2;


global OptiInputStruct

blo=OptiInputStruct.blo;
% overwrite parameters with blowing configuration

blo.active=true;                 %  activate blowing
blo.L= {[z(2)]*OptiInputStruct.prf.c;              %  length of blowing area: suction side
        [0.1]*OptiInputStruct.prf.c;};            %  length of blowing area: pressure side
blo.x= {xm*OptiInputStruct.prf.c;          %  midpoint of blowing area
        0.1*OptiInputStruct.prf.c;};             
blo.A= {[z(3)]*OptiInputStruct.flo.Uinfty;
        [0.0]*OptiInputStruct.flo.Uinfty};
blo.pressureCor=false;        % include correction term for pressure

try
[sol,~,~,~,~,~]=airfoil(OptiInputStruct.prf,OptiInputStruct.flo,...
                        OptiInputStruct.tri,blo,OptiInputStruct.eng,OptiInputStruct);

%disp(['Upper ',num2str(xB(1)),' lower ',num2str(xB(2)), ' ratio ',num2str(sol.CL/sol.Cdrag)])
ratio=-sol.CL/sol.Cdrag;

%make sure unconvereged solution is not found as solution
if sol.residual>1e-3;
   ratio=0; 
end
catch % if there is any error
   ratio=0; 
end

end

%PlotProfile(prf,flo.wake,sol, 3);

