function [ lsg1 , ShapeParams ] = FalknerSkanInt( beta, reverseRegion,inverseFlow,doPlot)
% FALKNERSKANINT solves the Falkner-Skan equation using a shooting method with a Runge-Kutta integration method (4. order).
%                Integration for beta parameters is done using midpoint rule approximation
%                input: beta          - angle
%                       reverseRegion - if true solution with negative f''(0) for -0.199< beta < 0 
%                                           -> default set to false   
%                       InverseFlow   - if true flow direction is reversed (alpha=-1) 
%                                           -> default set to false    
%                       doPlot        - plots the solution
%                                           -> default set to false             

N=200;
grading=1.006;
itmax=30;


if nargin==1
    reverseRegion=false;
    doPlot=false;
    inverseFlow=false;
elseif nargin==2
    inverseFlow=false;
    doPlot=false;
elseif nargin==3
    doPlot=false;
end

% check which branch of beta occurs
if beta>=0
    etaEnd= 5;
    deltas=0.005;
elseif -0.199<= beta && beta<0
    etaEnd= 8;
    deltas=0.001;
else % beta < -0.199  -> seperation occurs, reverse flow
    N=320;
    itmax=70;
    etaEnd= 20;
    deltas=0.008;
end


ind= 1:N; % vector with node indices
dim=size(ind);

% get stepsize vector and diskretisiced coordinate vector
if grading==1 % no grading
    Deta_in=etaEnd/N;
    Deta=Deta_in*ones(dim);
else
    gn=grading.^ind;
    mgn=sum(gn);
    Deta_in= etaEnd/mgn;
    Deta= Deta_in*gn;
end
eta= [0, cumsum(Deta)];


% table values for b -> interpolate to use as starting guess for s
% values from Schlichting, Cebecci and Asaithambi
bet=[-2     ,-1.6    ,-1.2   ,-1      ,-0.6    ,-0.199,-0.1988 ,-0.18    , -0.15  ,-0.12   ,-0.1     , -0.05    , 0       , 0.1     , 0.2     , 0.3   , 0.4    , 0.5    , 0.6      , 0.7   , 0.8   , 0.9   , 1      ,1.2     , 1.6    ,20/11  ,200/101, 2];
st =[-1.5823,-1.40457,-1.2015,-1.08638,-0.81202,0     ,0.005216, 0.128636,0.216362,0.281761,0.3192733, 0.4003238,0.4696005, 0.587037, 0.686711, 0.7681, 0.85442, 0.92768, 0.9958366, 1.0590, 1.1202, 1.1777, 1.23261,1.335724,1.521516,1.61399,1.67940, 1.687218];
% beta between -0.199 and zero solution with reverse flow and solution without reverse flow possible
if reverseRegion
    ind= find(bet>-0.199 & bet<0);
    btt=[-0.18    , -0.15    , -0.12    , -0.1    ,-0.05 ,-0.025];
    stt=[-0.097692, -0.133421, -0.142935,-0.140546,-0.108,-0.074];
    bet(ind)=[];  bet= [bet(1:ind(1)-1),btt,bet(ind(1):end)];
    st(ind)=[];  st= [st(1:ind(1)-1) ,stt,st(ind(1):end)];
end

% figure
% plot(bet(2:end),st(2:end));
% ylabel('$ \partial_{\eta \eta} f(0) $') 
% xlabel('$ \beta $') 
% title('interpolated values s(\beta) from table -> used as initial guess')

in=find(bet>=beta,1);
if in==1 
    s_approx=st(in);
else   
    s_approx= (beta-bet(in-1))/(bet(in)-bet(in-1))* st(in) + (bet(in)-beta)/(bet(in)-bet(in-1))* st(in-1);
end

% initial guess for value of f''(0)
s1=s_approx;

Initial= [0; 0; s1];
% integrate initial value Problem with Runge-Kutta method
[lsg1, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow);


% if solution diverges try other starting values until no divergence occurs
%-----------------------------------------------------------------
k=1;direction=1;
while diverged && k<itmax
    stmp=s1 + k*deltas;
    Initial= [0; 0; stmp];
    [lsg1, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow);
    if diverged
        stmp=s1 - k*deltas;
        Initial= [0; 0; stmp];
        [lsg1, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow);
    end
    k=k+1;
end
if diverged
   disp('no convergence at max iteration');
   return;
elseif k>1
    %disp(['first non diverging solution s=',num2str(stmp), ' it=', num2str(k-1)]);  
    del= s1-stmp;
    if del<-1e-4; direction=-1; else direction=1; end
    s1=stmp;
end

% f'(infty) -> condition for right boundary value
dfinf1=lsg1(2,end);

%-----------------------------------------------------------------

% second solution to start secant method
s2=s1 + direction*deltas;
Initial= [0; 0; s2];
[lsg2, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow);
% if solution diverges try other starting values until no divergence occurs
k=1;
while diverged && k<itmax
    s2=s1 + direction*k*deltas;
    Initial= [0; 0; s2];
    [lsg2, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow);
    k=k+1;
end
dfinf2=lsg2(2,end);
% approximated derivate for Newton method
dg_ds = (dfinf2-dfinf1) /deltas;
g= dfinf2 -1;

deltas= -g/dg_ds;


%--------- shooting method -------------------------

intersect=false;
k=0; minM=-1; minP=1;
% shooting method with Newton method for staring value s
while abs(deltas)> 1e-7 && k<itmax
    s2=s1+deltas;
    Initial= [0; 0; s2];
    [lsg2, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow);
    
    if ~ diverged
        % f'(infty)
        dfinf2=lsg2(2,end);
        
        % approximated derivate for Newton method
        dg_ds = (dfinf2-dfinf1) /deltas;
        g= dfinf2 -1;
        
        % intersection area
        if g < 0 && g > minM
            minM=g;
            sM=s2; dM=dfinf2;
        elseif g > 0 && g < minP
            minP=g;
            sP=s2; dP=dfinf2;
        end

        deltas= -g/dg_ds;

        s1=s2;
        dfinf1=dfinf2;
        lsg1=lsg2;
    else
        if minM ~=-1 && minP ~=1; 
            intersect=true;
        end
        break;
    end
    k=k+1;
end

% intersection method
if intersect % Newton method diverging -> try intersection method
   s1=sM; s2=sP;
   v1=dM; v2=dP;
   k=0;
   while abs(del)>1e-3 && k< itmax
       stmp=0.5*(s1+s2);
       Initial= [0; 0; stmp];
       [lsg2, diverged]=RungeKuttaFSeq( beta, Initial, Deta ,inverseFlow); 
       dfinf2=lsg2(2,end);

       dd2= v2-dfinf2;       
       if  dd2 < 0;  
           s1=stmp; v1=dfinf2; 
       else
           s2=stmp; v2=dfinf2; 
       end

       k=k+1;
   end
elseif diverged
    disp(['solution diverged, f´(infty)=', num2str(dfinf1)])
end


% calculate all Parameters

df= lsg1(2,:);
ddf=lsg1(3,:);

% f'''(0)=-beta
% approximation with forward differences
d3fw=( lsg1(3,2)-lsg1(3,1) )/(Deta(1) );

% error estimate for angle beta= - f'''(0)
e_abs= (beta+d3fw);
e_rel= abs(e_abs/beta);


% displacement thickness factor delta_1=beta1*delta_N
beta1= eta(end)-lsg1(1,end);

% momentum thickness factor delta_2=beta2*delta_N
% -> beta2= (f''(0)-beta*beta1 ) /( beta+1)
if beta ~= -1
    beta2= ( s1 - beta*beta1 )/(beta+1);
else % calculate via integration if singular
    g2= ((1 - df).*df);
    beta2= sum( 0.5*(g2(1:end-1)+g2(2:end)) .*Deta );
end

%midpoint rule for integration
% kinetic energie thickness factor delta_3=beta3*deltaN -> beta3= int (1-f'^2)f'
g3= ((1 - df.^2).*df);
beta3= sum( 0.5*(g3(1:end-1)+g3(2:end)) .*Deta );

% % Dissipation factor
gD=ddf.^2;
betaD= sum( 0.5*(gD(1:end-1)+gD(2:end)) .*Deta );

H32=beta3/beta2;
H12=beta1/beta2;


% Cf *ReT = 2*f''(0)*beta2
Cf=2*s1*beta2;

% CD *ReT = 2*f''(0)*beta2
CD=beta2*betaD;

ShapeParams.s=s1;
ShapeParams.beta1=beta1;
ShapeParams.beta2=beta2;
ShapeParams.beta3=beta3;
ShapeParams.betaD=betaD;
ShapeParams.CD=CD;
ShapeParams.Cf=Cf;
ShapeParams.H32=H32;
ShapeParams.H12=H12;

% Plotting
if doPlot
    % % str={['initial value s=\partial_{\eta \eta} f(0)=' num2str(s1)]...
    % %     ,['error   e= \beta + \partial_{\eta \eta \eta}f(0)=' num2str(e_rel*100) '%']};
    % 

    str={ ['initial value s=\partial_{\eta \eta} f(0)=' num2str(s1)]...
        ,['error   e= \beta + \partial_{\eta \eta \eta}f(0)=' num2str(e_rel*100) '%'],...
        ['H_{12}=', num2str(H12)] ,['H_{32}=', num2str(H32)],['C_{f}Re_{\delta_2}=', num2str(Cf)],['C_{D}Re_{\delta_2}=', num2str(CD)]  ...
        ['\beta_1=',num2str(beta1)], ['\beta_2=',num2str(beta2)],  ['\beta_3=',num2str(beta3)], ['\beta_D=',num2str(betaD)] };

    if beta>-0.199; tst=34; else tst=60; end
    
    figure 
    hold on
    plot(eta,lsg1(1,:));
    plot(eta,lsg1(2,:));
    plot([eta(1), eta(end)],[1,1]);
    %text(etaEnd/10,lsg1(1,end-10),str);
    text(etaEnd/10,lsg1(1,end-tst),str);
    title(['Falkner-Skan solution for \beta=', num2str(beta)]);
    xlabel(' \eta ') ;
    legend('$f$' ,'$\partial_{\eta} f$' ,'location','northeast'); 
end

end


% % midpoint rule for integration
% df= 0.5*(lsg1(2,1:end-1) + lsg1(2,2:end));
% ddf= 0.5*(lsg1(3,1:end-1) + lsg1(3,2:end));
% 
% % kinetic energie thickness factor delta_3=beta3*deltaN -> beta3= int (1-f'^2)f'
% beta3t= sum( ((1 - df.^2).*df).*Deta );
% 
% % % Dissipation factor
% betaDt= sum( ddf.^2.*Deta );





% m= H^2* (0.058*(H-4)^2/(H-1) -0.068) / ( 6.54*H -14.07);
% dm_dH= 2*m/H + 0.116*H^2*(H-4)/( (H-1)*( 6.54*H -14.07) ) - 0.058*(H-4)^2/ ( (H-1)^2*(6.54*H -14.07) )...
%             -6.54*m/( 6.54*H -14.07);
% 
% beta= 2*m / (m+1);
% dbeta_dm= (2-beta)/(m+1);

% m= beta/(2-beta);
% dm_dbeta= (1+m)/(2-beta);

% RBen f(0)=0 f'(0)=0  f'(infty)=1

% Initial guess for f''(0) -> s



% % initial solution guess
% if H<=3
%     % Sinus curve
%     omega= pi/10;
%      
% else
% 
% % Zustandsraum:   z1:=f , z2:= f' , z3:=f''
% 
% diff_z1= z2;
% diff_z2= z3;
% diff_z3= -z1.*z3 - beta* (1-z2.^2);
% 
% rhs=zeros(3*N,1);
% rhs(1:3:end-2)=diff_z1;
% rhs(2:3:end-1)=diff_z2;
% rhs(3:3:end  )=diff_z3;
% 
% % Jacobimatrix
% 
% % J= [  0 ,    1     , 0;...
% %       0 ,    0     , 1;...
% %      -z3, 2*beta*z2, -z1];
% 
% 
% Jges= zeros(3*N,3*N);
% for i=1:N
%  Ji=[  0 ,    1     , 0;...
%        0 ,    0     , 1;...
%       -z3(i), 2*beta*z2(i), -z1(i)];
%   Jges( 3*i-2:3*i ,3*i-2:3*i)= Ji; 
% end
 
 




