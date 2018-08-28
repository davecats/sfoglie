function [sol,prf,flo,tri,blo,eng,InvRef]=airfoil(prf,flo,tri,blo,eng,InvRef)
% AIRFOIL   calculate viscous solution
%           InvRef: if Inviscid solution already exists from previous calculations
%                       -> commit InvRef struct with fields prf, flo, CoeffMatrix, Uinv and sges
%                       -> algorithm checks whether profile geometry or angle of attack changes,
%                          if this is not the case, the previous inv sol can be used to save calc time 

%----------------------------------------------------------------------------
%%  INVISCID SOLUTION
%----------------------------------------------------------------------------

if nargin==6 % InvRef is commited
    % check whether relevant parameters have changes or not
    Nochange= sum( prf.naca == InvRef.prf.naca )==4 && prf.M==InvRef.prf.M && prf.pmode ==InvRef.prf.pmode...
              &&  prf.c==InvRef.prf.c && flo.alfa ==InvRef.flo.alfa && flo.Uinfty ==InvRef.flo.Uinfty...
              &&  prf.noSkew==InvRef.prf.noSkew && prf.sharpTE==InvRef.prf.sharpTE; 
else
    Nochange=false;
end

if Nochange
    % Inviscid solution already exists, no need to recalculate it
    try
        % assign values
        prf=InvRef.prf;
        tmp=flo.nkrit;
        flo=InvRef.flo;
        flo.nkrit=tmp;
        CoeffMatrix=InvRef.CoeffMatrix;
        Uinv=InvRef.Uinv;
        sges=InvRef.sges;
        disp('Inviscid solution already avaivable -> only do viscous part')
    catch % if struct is not complete still do inv. solution process
        Nochange=false;
    end
end
    
if ~Nochange
    % calculates the inviscid solution with the boundary element method
    [prf,flo,CoeffMatrix,Uinv,sges ] = InviscidSolution(prf,flo,eng);
    disp('Inviscid solution was calculated -> start BL calcs')
     % assign InvRef for output to be used for other simulations
     InvRef.prf=prf;
     InvRef.flo=flo;
     InvRef.CoeffMatrix=CoeffMatrix;
     InvRef.Uinv=Uinv;
     InvRef.sges=sges;
end


% if only inviscid solution is required
if flo.invisc
    sol.invisc=true;
    sol.U=Uinv;
    sol.Cp= 1-(Uinv/flo.Uinfty).^2;
    sol.CL= getCL(prf,Uinv,flo.alfa,flo.Uinfty);
    return
end


%----------------------------------------------------------------------------
%%  VISCOUS SOLUTION
%----------------------------------------------------------------------------

% blowing velocity vector
blo.Vb=zeros(size(Uinv));
% get nodes with blowing and set the blowing velocity vector
blo=addBlowing(blo,prf);
% initial solution guess
sol = GetInitialSolution( prf,flo, tri, eng, Uinv, blo.Vb, 2);
%  coupled boundary layer and potential flow solution
[sol, prf]=NewtonEq( prf,flo,eng,sol,CoeffMatrix.D,Uinv,blo.pressureCor);

end

%PlotProfile(prf,flo.wake,sol, 3);







