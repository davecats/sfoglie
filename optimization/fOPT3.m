function  ratio  = fOPT3(z)
%FOPT3 function for optimization two regions over the arclength.
%      Parameters: s1_start, s1_end, s2_start, s2_end, Vb1, Vb2


% the algorithm can not deal with variable constraints so that kind of dirty
% trick has to be used. setting ratio to 0 will make the algorithm go away
% from solutions where blowing end point would be out of the bounds
if z(2)< z(1) || z(3)< z(2) || z(4)< z(3)
    ratio=0; 
    return;
end


global OptiInputStruct

blo=OptiInputStruct.blo;
% overwrite parameters with blowing configuration
blo.pressureCor=false;        % include correction term for pressure

blo.active=true;
blo.ArcLengthMode=true;

% arclengths in percent of smax
blo.s_change=        [z(1),z(2),z(3),z(4)]; 
% intensitie for s_i< s < s_i+1
blo.NewA= OptiInputStruct.flo.Uinfty*[z(5),0,z(6),0];



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

