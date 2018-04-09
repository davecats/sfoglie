%  set parameters and runs the calculation script


% general setup

p.NACA=[4 4 1 2];
p.Re= 1e5;
p.alfa= 0; % angle of attack in degrees
p.Uinfty=1;

% profile
%-------------------------
p.M=90; % diskretication points on x-axis -> number of nodes N= 2M - 1
p.NoSkew  =false;        % if true profile skewness neclected
p.sharpTE =false;        % if true sharp Trailing edge   
p.c=1; %chordlength



% Transition
%-------------------------
p.nkrit=9;

p.trip=[ false;...       % tripping on suction side 
         false];         % tripping on pressure side 

%tripping location
p.xtrip=[ 0.04;...
          1.0 ]*p.c;

      
% Blowing
%-------------------------     
p.withBlowing=[false;...  % blowing on suction side  
               false];   % blowing on pressure side 
           
% blowing region
% startpoint    
p.xBstart= [0.25;...
            0.25]* p.c;
% end point     
p.xBend  = [0.5;...
            0.5]* p.c;
      
% blowing intensity      
p.intensity=[0.001;...
             0.001]* p.Uinfty;

p.pressureCor=false;% true;%

p.it=24;

Calc



return;

% loop over alfas + write out   for CL over alfa plots
for k=0:12

    p.alfa= k;
    disp(['ALFA=',num2str(p.alfa)])
    disp('--------------------------------------------')
    
    
    % evaluation
    Calc

    if k==0
        CP=sol.Cp;
        TAU=sol.tau;
        CL   =sol.CL ;
        Cnu  =sol.Cnu;
        Cdrag=sol.Cdrag;
        alph=p.alfa;
    else
        CP=[CP,sol.Cp];
        TAU=[TAU,sol.tau];
        CL=[CL,sol.CL];
        Cnu=[Cnu,sol.Cnu];
        Cdrag=[Cdrag,sol.Cdrag];
        alph=[alph,p.alfa];
    end

end


 tmp=Re; i=0;
while tmp>1
    tmp=tmp/10;
    i=i+1;  
end
first=round(tmp*10);
str1=['Re',num2str(first),'e',num2str(i-1)];

if trip(1)
    str2='_Trip';
else
    if nkrit>1
        str2=['_N',num2str(round(nkrit))];
    else
        tmp=100*nkrit;
        stmp=num2str(tmp);
        if strcmp(stmp(2),'0'); stmp=stmp(1); end
        str2=['_N0',stmp];
    end

end
str=[str1,str2];

if p.withBlowing(1) || p.withBlowing(1)
   str=['Blow_',str] ;
end

glob=[alph;CL;Cdrag;Cnu];
dlmwrite([str,'_cp'],CP)
dlmwrite([str,'_tau'],TAU)
dlmwrite([str,'_glob'],glob)





