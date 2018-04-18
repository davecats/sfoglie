function [ J, rhs] = JacobiM( prf,flo,blo,eng,sol,Unew,D)
%JACOBIM        calculates the Jacobimatrix and right hand side for the Newton system 
%               for solution vector z=[T1,..,TN, C1,.., CN,m1,..,mN]^T 


nu=flo.nu;
dim=size(Unew);
Nges=prf.N + flo.wake.N;

% panel length vector
h=[prf.panels.L,flo.wake.L]';


%prescribe
f1=zeros(dim);
f2=zeros(dim);
f3=zeros(dim);
df1_dT=zeros(Nges,2);
df1_dD=zeros(Nges,2);
df1_dU=zeros(Nges,2);
%df1_dC=zeros(Nges,2);
df1_ds=zeros(Nges,2);
df2_dT=zeros(Nges,2);
df2_dD=zeros(Nges,2);
df2_dU=zeros(Nges,2);
df2_dC=zeros(Nges,2);
df2_ds=zeros(Nges,2);
df3_dT=zeros(Nges,2);
df3_dD=zeros(Nges,2);
df3_dU=zeros(Nges,2);
df3_dC=zeros(Nges,2);
df3_ds=zeros(Nges,2);
%                                        laminar part of Boundary
%------------------------------------------------------------------------------------------------------------------------

% suction side
%------------------------------------------------------------

ind=(prf.Nle-1 :-1: sol.iTran(1)+1); % laminar node indizes starting from LE
indM=ind(2:end); % panel indizes
s1=prf.sU(ind(1:end-1))';


if length(ind)>1 
    [ f1l,f2l,der] =JacobiLam(sol.T(ind),sol.U(ind),sol.Vb(ind),sol.HK(ind),sol.Ret(ind),h(indM),s1,nu);


    f1(indM)=f1l;
    f2(indM)=f2l;
    df1_dT(indM,:)=der.df1_dT;
    df1_dD(indM,:)=der.df1_dD;
    df1_dU(indM,:)=der.df1_dU;
    df1_ds(indM,:)=der.df1_ds;
    df2_dT(indM,:)=der.df2_dT;
    df2_dD(indM,:)=der.df2_dD;
    df2_dU(indM,:)=der.df2_dU;
    df2_ds(indM,:)=der.df2_ds;


    % Amplification Equation
    n=sol.c(ind); 

     [ f3(indM),df3_dC(indM,:), df3_dT(indM,:),df3_dD(indM,:),df3_dU(indM,:),df3_ds(indM,:) ] =...
                AmplificationEquation(flo,n,sol.T(ind),sol.U(ind),sol.HK(ind),sol.Ret(ind),h(indM)  );
end

% pressure side
%------------------------------------------------------------

ind=(prf.Nle: sol.iTran(2)-1);
indM=ind(1:end-1);
s1=prf.sL(indM-prf.Nle+1)';

if length(ind)>1 
    [ f1l,f2l,der] = JacobiLam(sol.T(ind),sol.U(ind),sol.Vb(ind),sol.HK(ind),sol.Ret(ind),h(indM),s1,nu);

    f1(indM)=f1l;
    f2(indM)=f2l;
    df1_dT(indM,:)=der.df1_dT;
    df1_dD(indM,:)=der.df1_dD;
    df1_dU(indM,:)=der.df1_dU;
    df1_ds(indM,:)=der.df1_ds;
    df2_dT(indM,:)=der.df2_dT;
    df2_dD(indM,:)=der.df2_dD;
    df2_dU(indM,:)=der.df2_dU;
    df2_ds(indM,:)=der.df2_ds;

    % Amplification Equation
    n=sol.c(ind); 
     [ f3(indM),df3_dC(indM,:), df3_dT(indM,:),df3_dD(indM,:),df3_dU(indM,:),df3_ds(indM,:) ] =...
                AmplificationEquation(flo,n,sol.T(ind),sol.U(ind),sol.HK(ind),sol.Ret(ind),h(indM)  );
end

%                               turbulent part of Boundary
%------------------------------------------------------------------------------------------------------------------------

% suction side
%------------------------------------------------------------
ind=(sol.iTran(1) :-1: 1); %turbulent node indizes
if length(ind)>1 % check for turbulent flow
    indM=ind(2:end); % indices for midpoint values -> without transition panel
    s1=prf.sU(ind(1:end-1))';
    
    [ f1t,f2t,f3t,der ] = JacobiTurb(sol.D(ind),sol.T(ind),sol.c(ind),sol.U(ind),sol.Vb(ind),sol.HK(ind),sol.Ret(ind),h(indM),s1,false,false,nu);

    f1(indM)=f1t;
    f2(indM)=f2t;
    f3(indM)=f3t;

    df1_dT(indM,:)=der.df1_dT;
    df1_dD(indM,:)=der.df1_dD;
    df1_dU(indM,:)=der.df1_dU;
    df1_ds(indM,:)=der.df1_ds;
    df2_dT(indM,:)=der.df2_dT;
    df2_dD(indM,:)=der.df2_dD;
    df2_dU(indM,:)=der.df2_dU;
    df2_dC(indM,:)=der.df2_dC;
    df2_ds(indM,:)=der.df2_ds;
    df3_dT(indM,:)=der.df3_dT;
    df3_dD(indM,:)=der.df3_dD;
    df3_dU(indM,:)=der.df3_dU;
    df3_dC(indM,:)=der.df3_dC;
    df3_ds(indM,:)=der.df3_ds;
end


% pressure side
%------------------------------------------------------------
ind=(sol.iTran(2):prf.N); %turbulent node indizes
if length(ind)>1
    indM=ind(1:end-1); % indices for midpoint values -> without transition panel
    s1=prf.s(sol.iTran(2):end-1)';
    
    [ f1t,f2t,f3t,der ] = JacobiTurb(sol.D(ind),sol.T(ind),sol.c(ind),sol.U(ind),sol.Vb(ind),sol.HK(ind),sol.Ret(ind),h(indM),s1,false,false,nu);
    
    f1(indM)=f1t;
    f2(indM)=f2t;
    f3(indM)=f3t;

    df1_dT(indM,:)=der.df1_dT;
    df1_dD(indM,:)=der.df1_dD;
    df1_dU(indM,:)=der.df1_dU;
    df1_ds(indM,:)=der.df1_ds;
    df2_dT(indM,:)=der.df2_dT;
    df2_dD(indM,:)=der.df2_dD;
    df2_dU(indM,:)=der.df2_dU;
    df2_dC(indM,:)=der.df2_dC;
    df2_ds(indM,:)=der.df2_ds;
    df3_dT(indM,:)=der.df3_dT;
    df3_dD(indM,:)=der.df3_dD;
    df3_dU(indM,:)=der.df3_dU;
    df3_dC(indM,:)=der.df3_dC;
    df3_ds(indM,:)=der.df3_ds;
end

% flo.wake
%------------------------------------------------------------
ind=(prf.N+1:prf.N+flo.wake.N);
indM=ind(1:end-1);
s1=flo.wake.s(1:end-1)'+prf.sL(end);


[ f1t,f2t,f3t,der ] = JacobiTurb(sol.D(ind),sol.T(ind),sol.c(ind),sol.U(ind),sol.Vb(ind),sol.HK(ind),sol.Ret(ind),h(indM),s1,true,flo.wake.gap,nu);

f1(indM)=f1t;
f2(indM)=f2t;
f3(indM)=f3t;

df1_dT(indM,:)=der.df1_dT;
df1_dD(indM,:)=der.df1_dD;
df1_dU(indM,:)=der.df1_dU;
df2_dT(indM,:)=der.df2_dT;
df2_dD(indM,:)=der.df2_dD;
df2_dU(indM,:)=der.df2_dU;
df2_dC(indM,:)=der.df2_dC;
df2_ds(indM,:)=der.df2_ds;
df3_dT(indM,:)=der.df3_dT;
df3_dD(indM,:)=der.df3_dD;
df3_dU(indM,:)=der.df3_dU;
df3_dC(indM,:)=der.df3_dC;
df3_ds(indM,:)=der.df3_ds;


%                               Transition panel
%------------------------------------------------------------------------------------------------------------------------

% suction side
%------------------------------------------------------------

ind= (sol.iTran(1)+1:-1:sol.iTran(1));
indM=sol.iTran(1);
n= [sol.c(ind(1));sol.tran.n2(1)];

s1=prf.sU(sol.iTran(1)+1);
% hl=sol.tran.Llam(1);
% st=sl+hl;

if ~sol.Tripping(1) % free Transition
    [ fT, der,~ ] = TransitionEQ( flo, eng, n, sol.T(ind), sol.D(ind),sol.U(ind),sol.Vb(ind),s1,h(indM), sol.c(ind(2)) );   
else % forced Transition 
    [ fT, der,~ ] = TransitionEQ( flo, eng, n, sol.T(ind), sol.D(ind),sol.U(ind),sol.Vb(ind),s1,h(indM), sol.c(ind(2)) ,sol.sT(1),true);   
    % add forced transition derivate
    der.df1_ds(2)=der.df1_ds(2) + der.df1_dsF;
    der.df2_ds(2)=der.df2_ds(2) + der.df2_dsF;
    der.df3_ds(2)=der.df3_ds(2) + der.df3_dsF;
end 




f1(indM)=fT(1);
f2(indM)=fT(2);
f3(indM)=fT(3);

df1_dT(indM,:)=der.df1_dT;
df1_dD(indM,:)=der.df1_dD;
df1_dU(indM,:)=der.df1_dU;
df1_ds(indM,:)=der.df1_ds;
df2_dT(indM,:)=der.df2_dT;
df2_dD(indM,:)=der.df2_dD;
df2_dU(indM,:)=der.df2_dU;
df2_dC(indM,:)=der.df2_dC;
df2_ds(indM,:)=der.df2_ds;
df3_dT(indM,:)=der.df3_dT;
df3_dD(indM,:)=der.df3_dD;
df3_dU(indM,:)=der.df3_dU;
df3_dC(indM,:)=der.df3_dC;
df3_ds(indM,:)=der.df3_ds;



% pressure side
%------------------------------------------------------------

ind= (sol.iTran(2)-1:sol.iTran(2));
indM=sol.iTran(2)-1;
n= [sol.c(ind(1));sol.tran.n2(2)];

s1=prf.sL(sol.iTran(2)-prf.Nle);
% hl=sol.tran.Llam(2);
% st=sl+hl;

if ~sol.Tripping(2) % free Transition
    [ fT, der,~ ] = TransitionEQ( flo, eng, n, sol.T(ind), sol.D(ind),sol.U(ind),sol.Vb(ind),s1,h(indM), sol.c(ind(2)) );  
else % forced Transition
    [ fT, der,~ ] = TransitionEQ( flo, eng, n, sol.T(ind), sol.D(ind),sol.U(ind),sol.Vb(ind),s1,h(indM), sol.c(ind(2)) ,sol.sT(2),true);  
    % add forced transition derivate  
    der.df1_ds(2)=der.df1_ds(2) + der.df1_dsF;
    der.df2_ds(2)=der.df2_ds(2) + der.df2_dsF;
    der.df3_ds(2)=der.df3_ds(2) + der.df3_dsF;  
end  



f1(indM)=fT(1);
f2(indM)=fT(2);
f3(indM)=fT(3);

df1_dT(indM,:)=der.df1_dT;
df1_dD(indM,:)=der.df1_dD;
df1_dU(indM,:)=der.df1_dU;
df1_ds(indM,:)=der.df1_ds;
df2_dT(indM,:)=der.df2_dT;
df2_dD(indM,:)=der.df2_dD;
df2_dU(indM,:)=der.df2_dU;
df2_dC(indM,:)=der.df2_dC;
df2_ds(indM,:)=der.df2_ds;
df3_dT(indM,:)=der.df3_dT;
df3_dD(indM,:)=der.df3_dD;
df3_dU(indM,:)=der.df3_dU;
df3_dC(indM,:)=der.df3_dC;
df3_ds(indM,:)=der.df3_ds;

%------------------------------------------------------------------------------------------------------------------------

% switch to m as variable
%----------------------------------------------------

% save the Indizes of equations and the corresponding values of 1 and 2 nodes

% node indizes
% suction side -> 1 and 2 are "switched" when walking in global arc length direction
Node1= (2:prf.Nle-1);
Node2= (1:prf.Nle-2);
% pressure side
Node1=[Node1, (prf.Nle:prf.N-1)];
Node2=[Node2, (prf.Nle+1:prf.N)];
% flo.wake
Node1=[Node1, (prf.N+1:prf.N+flo.wake.N-1)];
Node2=[Node2, (prf.N+2:prf.N+flo.wake.N)  ];

% Equation index -> no FDM equation for LE panel, TE panel and last flo.wake node -> 9 additional Equations neccessary
EQ = [(1:prf.Nle-2),(prf.Nle:prf.N-1),(prf.N+1:prf.N + flo.wake.N-1) ];

% different order for Equation System
EQsys= [1:prf.Nle-2, prf.Nle+1:prf.N, prf.N+2:prf.N + flo.wake.N];

% total derivates of equations in respect to m   

% Contribution of U
dU1_dm = D(Node1,:); 
dU2_dm = D(Node2,:); 

% Contribution of DI
tmp=ones(1,length(sol.U));
dD_dm=diag( 1./sol.U ) - ( sol.D./(sol.U)*tmp ).* D;


dD1_dm=dD_dm(Node1,:);
dD2_dm=dD_dm(Node2,:);

% stagnation point change
dss_dU1=prf.dsLE_dG1*ones(Nges,1);dss_dU1(prf.Nle:end)=-dss_dU1(prf.Nle:end);% switch sign for pressure side and flo.wake
dss_dU2=prf.dsLE_dG2*ones(Nges,1);dss_dU2(1:prf.Nle-1)=-dss_dU2(1:prf.Nle-1);% switch sign for suction side


% total derivates in respect to m
df1_dm=   (df1_dD(EQ,1)*tmp) .*dD1_dm + (df1_dD(EQ,2)*tmp) .*dD2_dm ...
        + (df1_dU(EQ,1)*tmp) .*dU1_dm + (df1_dU(EQ,2)*tmp) .*dU2_dm ...
        + ( df1_ds(EQ,1) + df1_ds(EQ,2) )*ones(size(D(1,:))) ...
           .* (dss_dU1(EQ)*D(prf.Nle-1,:) + dss_dU2(EQ)*D(prf.Nle,:) );


df2_dm=   (df2_dD(EQ,1)*tmp) .*dD1_dm + (df2_dD(EQ,2)*tmp) .*dD2_dm ...
        + (df2_dU(EQ,1)*tmp) .*dU1_dm + (df2_dU(EQ,2)*tmp) .*dU2_dm ...      
        + ( df2_ds(EQ,1) + df2_ds(EQ,2) )*ones(size(D(1,:))) ...
            .* (dss_dU1(EQ)*D(prf.Nle-1,:) + dss_dU2(EQ)*D(prf.Nle,:) );
    
    
df3_dm=   (df3_dD(EQ,1)*tmp) .*dD1_dm + (df3_dD(EQ,2)*tmp) .*dD2_dm ...
        + (df3_dU(EQ,1)*tmp) .*dU1_dm + (df3_dU(EQ,2)*tmp) .*dU2_dm ...
        + ( df3_ds(EQ,1) + df3_ds(EQ,2) )*ones(size(D(1,:))) ...
            .* (dss_dU1(EQ)*D(prf.Nle-1,:) + dss_dU2(EQ)*D(prf.Nle,:) );
      
 
% T-part 
% EQ 1
JT1=zeros(Nges,Nges);
JT1(EQsys,Node1)=diag(df1_dT(EQ,1));   
JT1(EQsys,Node2)=JT1(EQsys,Node2) + diag(df1_dT(EQ,2));  

% EQ 2
JT2=zeros(Nges,Nges);
JT2(EQsys,Node1)=diag(df2_dT(EQ,1));   
JT2(EQsys,Node2)=JT2(EQsys,Node2) + diag(df2_dT(EQ,2));     
% EQ 3
JT3=zeros(Nges,Nges);
JT3(EQsys,Node1)=diag(df3_dT(EQ,1));   
JT3(EQsys,Node2)=JT3(EQsys,Node2) + diag(df3_dT(EQ,2));   

% C-part
% EQ 2
JC2=zeros(Nges,Nges);
JC2(EQsys,Node1)=diag(df2_dC(EQ,1));   
JC2(EQsys,Node2)=JC2(EQsys,Node2) + diag(df2_dC(EQ,2));     
% EQ 3
JC3=zeros(Nges,Nges);
JC3(EQsys,Node1)=diag(df3_dC(EQ,1));   
JC3(EQsys,Node2)=JC3(EQsys,Node2) + diag(df3_dC(EQ,2));   


% total matrix
J= zeros(3*Nges,3*Nges);
EQ1= 3*EQsys - 2*ones(size(EQ));
EQ2= 3*EQsys - ones(size(EQ));
EQ3= 3*EQsys ;

% Model for correction of the pressure Term   
%------------------------------------------------------------------------------------------

 Blow=find(sol.Vb~=0);
 if ~isempty(Blow) && blo.pressureCor
     BlowU= Blow(Blow<prf.Nle-1);
     BlowL= Blow(Blow>prf.Nle);
     
     if ~isempty(BlowU)
        sU_b=prf.sU(BlowU);
        PrU = PressureCorrect(sU_b(end),sU_b(1),prf.sU(end:-1:1)',sol.Vb(prf.Nle-1:-1:1),sol.U(prf.Nle-1:-1:1));
        PrU=PrU(end:-1:1);
        Upper= PrU./sol.T(1:prf.Nle-1);
        Upper= 0.5*(Upper(2:end)+Upper(1:end-1)).*prf.panels.L(1:prf.Nle-2)';
        
        f1(1:prf.Nle-2)=f1(1:prf.Nle-2) + Upper;
     end
     if ~isempty(BlowL)
        sL_b=prf.sL(BlowL-prf.Nle+1);
        PrL = PressureCorrect(sL_b(1),sL_b(end),prf.sL',sol.Vb(prf.Nle:prf.N),sol.U(prf.Nle:prf.N));
        Lower=PrL./sol.T(prf.Nle:prf.N);
        Lower= 0.5*(Lower(2:end)+Lower(1:end-1)).*prf.panels.L(prf.Nle:prf.N-1)';
    
        f1(prf.Nle:prf.N-1)=f1(prf.Nle:prf.N-1) + Lower;
     end
 end
%------------------------------------------------------------------------------------------

f=zeros(3*Nges,1);
f(EQ1)=f1(EQ); f(EQ2)=f2(EQ); f(EQ3)=f3(EQ);

J(EQ1,1:Nges)= JT1(EQsys,:);
J(EQ2,1:Nges)= JT2(EQsys,:);
J(EQ3,1:Nges)= JT3(EQsys,:);
J(EQ2,Nges+1:2*Nges)= JC2(EQsys,:);
J(EQ3,Nges+1:2*Nges)= JC3(EQsys,:);
J(EQ1,2*Nges+1:3*Nges)= df1_dm;
J(EQ2,2*Nges+1:3*Nges)= df2_dm;
J(EQ3,2*Nges+1:3*Nges)= df3_dm;


%-------------- right Hand side  -------------------------------------------------------
% changes forced by new velocity Unew= Uold + D*dm
deltaU=  sol.U - Unew ;    % negativ changes in U 
DDS   = -sol.D./sol.U .*deltaU; % negativ changes in DI


% incorporate forced changes for velocity and displacement thickness
deltaU1= deltaU(Node1);
deltaU2= deltaU(Node2);
DDS1=DDS(Node1);
DDS2=DDS(Node2);

add=zeros(size(f));

% incorporate forced Leading edge arc length change
stagn= dss_dU1.*deltaU(prf.Nle-1) + dss_dU2.*deltaU(prf.Nle);


add(EQ1)=   df1_dU(EQ,1).*deltaU1 + df1_dU(EQ,2).*deltaU2 ...
          + df1_dD(EQ,1).*DDS1    + df1_dD(EQ,2).*DDS2 ...
          + (df1_ds(EQ,1) + df1_ds(EQ,2) ).*stagn(EQ) ;
add(EQ2)=   df2_dU(EQ,1).*deltaU1 + df2_dU(EQ,2).*deltaU2 ...
          + df2_dD(EQ,1).*DDS1    + df2_dD(EQ,2).*DDS2 ...
          + (df2_ds(EQ,1) + df2_ds(EQ,2) ).*stagn(EQ) ;
add(EQ3)=   df3_dU(EQ,1).*deltaU1 + df3_dU(EQ,2).*deltaU2 ...
          + df3_dD(EQ,1).*DDS1    + df3_dD(EQ,2).*DDS2 ...
          + (df3_ds(EQ,1) + df3_ds(EQ,2) ).*stagn(EQ) ;

     
       
%--------------------------------------------------------------------------

% Equations to close System
% Index with empty equation
EMP= [3*(prf.Nle-1)-2,3*(prf.Nle-1)-1,3*(prf.Nle-1),3*prf.Nle-2,3*prf.Nle-1,3*prf.Nle, 3*(prf.N+1)-2,3*(prf.N+1)-1,3*(prf.N+1)];


% first suction side node
%---------------------------------------
[fs1,fs2, J1 ] = InitialNodeSys(sol.T(prf.Nle-1),sol.D(prf.Nle-1),sol.U(prf.Nle-1),sol.Vb(prf.Nle-1),prf.LE1,nu);
dfs1_dT= J1(1,1);
dfs2_dT= J1(2,1);
dfs1_dm= J1(1,2)*dD_dm(prf.Nle-1,:) + J1(1,3)*D(prf.Nle-1,:) ...
            + J1(1,4)*(dss_dU1(prf.Nle-1)*D(prf.Nle-1,:) + dss_dU2(prf.Nle-1)*D(prf.Nle,:) );
        
dfs2_dm= J1(2,2)*dD_dm(prf.Nle-1,:) + J1(2,3)*D(prf.Nle-1,:) ...
            + J1(2,4)*(dss_dU1(prf.Nle-1)*D(prf.Nle-1,:) + dss_dU2(prf.Nle-1)*D(prf.Nle,:) );
        
J(EMP(1), prf.Nle-1)=dfs1_dT; 
J(EMP(1), 2*Nges+1:end)=dfs1_dm; 
f(EMP(1))=fs1;
add(EMP(1))=   J1(1,3)*deltaU(prf.Nle-1) + J1(1,2)*DDS(prf.Nle-1) + J1(1,4) .*stagn(prf.Nle-1) ;


J(EMP(2), prf.Nle-1)=dfs2_dT; 
J(EMP(2), 2*Nges+1:end)=dfs2_dm; 
f(EMP(2))=fs2;
add(EMP(2))=   J1(2,3)*deltaU(prf.Nle-1) + J1(2,2)*DDS(prf.Nle-1) + J1(2,4).*stagn(prf.Nle-1) ;

% no Amplification change at initial points
J(EMP(3), Nges + prf.Nle-1)= 1;


% first pressure side node
%---------------------------------------
[fp1,fp2, J2 ] = InitialNodeSys(sol.T(prf.Nle),sol.D(prf.Nle),sol.U(prf.Nle),sol.Vb(prf.Nle),prf.LE2,nu);
dfp1_dT= J2(1,1);
dfp2_dT= J2(2,1);
dfp1_dm= J2(1,2)*dD_dm(prf.Nle,:) + J2(1,3)*D(prf.Nle,:)...
            + J2(1,4)*(dss_dU1(prf.Nle)*D(prf.Nle-1,:) + dss_dU2(prf.Nle)*D(prf.Nle,:) );
dfp2_dm= J2(2,2)*dD_dm(prf.Nle,:) + J2(2,3)*D(prf.Nle,:)...
            + J2(2,4)*(dss_dU1(prf.Nle)*D(prf.Nle-1,:) + dss_dU2(prf.Nle)*D(prf.Nle,:) );
        
J(EMP(4), prf.Nle)=dfp1_dT; 
J(EMP(4), 2*Nges+1:end)=dfp1_dm; 
f(EMP(4))=fp1;
add(EMP(4))=   J2(1,3)*deltaU(prf.Nle) + J2(1,2)*DDS(prf.Nle) + J2(1,4) .*stagn(prf.Nle) ;

J(EMP(5), prf.Nle)=dfp2_dT; 
J(EMP(5), 2*Nges+1:end)=dfp2_dm; 
f(EMP(5))=fp2;
add(EMP(5))=   J2(2,3)*deltaU(prf.Nle) + J2(2,2)*DDS(prf.Nle) + J2(2,4).*stagn(prf.Nle) ;


J(EMP(6), Nges + prf.Nle  )= 1;


% conditions for first flo.wake node
%---------------------------------------

% T_N+1=T_1 + T_N
J(EMP(7),1)      = 1;
J(EMP(7),prf.N)  = 1;
J(EMP(7),prf.N+1)=-1;


% f= D_N+1 - D_1 - D_N - L_TE = 0
fw= sol.D(prf.N+1) - sol.D(1) - sol.D(prf.N) - prf.gap;
dfw_dm= dD_dm(prf.N+1,:) -  dD_dm(1,:) - dD_dm(prf.N,:);

DDSw   = DDS(prf.N+1) - DDS(1) - DDS(prf.N);

J(EMP(8),2*Nges+1:end)=dfw_dm; 
f(EMP(8))= fw ;
add(EMP(8))=DDSw;


% f= - C_N+1 + (T_1*C_1 + T_N*C_N )/ (T_1+T_N) =0
CN1= ( sol.c(1)*sol.T(1) + sol.c(prf.N)*sol.T(prf.N) )/( sol.T(1)+sol.T(prf.N) );
dCN1_dC1= sol.T(1)    /( sol.T(1)+sol.T(prf.N) );
dCN1_dCN= sol.T(prf.N)/( sol.T(1)+sol.T(prf.N) );
dCN1_dT1= (sol.c(1)     -CN1)  /( sol.T(1)+sol.T(prf.N) );
dCN1_dTN= (sol.c(prf.N) -CN1)  /( sol.T(1)+sol.T(prf.N) );
J(EMP(9),1)=dCN1_dT1;
J(EMP(9),prf.N)=dCN1_dTN;
J(EMP(9),Nges+1)=dCN1_dC1;
J(EMP(9),Nges+prf.N)=dCN1_dCN;
J(EMP(9),Nges+prf.N+1)=-1; 

f(EMP(9))=-sol.c(prf.N+1) + CN1;


% combine right hand side
rhs=-f;
rhs=rhs + add;



end

