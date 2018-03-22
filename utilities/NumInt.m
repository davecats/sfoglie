function [ I,fint ] = NumInt( f,x, mode, cum )
%NumInt Does a numerical Integration where x are are nodes and f the corresponding values
%       mode=1   - uses Simpson law (default)
%       mode=2   - uses midpoint law
%       cum=true - returns the vector fint= int_a^x f dx


if nargin<3
    mode=1;
end
if nargin<4
    cum=false;
end

% make sure the dimensions fit
if size(x,1)<2;
   x=x'; 
end
if size(f,1)<2;
   f=f'; 
end

if mode==1
    %---------------------------------------------------------
    %           Simpson law
    %---------------------------------------------------------
     
     
    % check if number of intervalls is even
    res= mod(length(x)-1,2);
    
    if res==0 % even number simpson law is posible everywhere
        % nodes of each integration intervall
        x1=x(1:2:end-2); % left integration border
        x2=x(2:2:end-1); % midpoint node
        x3=x(3:2:end  ); % right integration border
        
        f1=f(1:2:end-2); % left integration border
        f2=f(2:2:end-1); % midpoint node
        f3=f(3:2:end  ); % right integration border
          
        % weigthening factors (integrals of Lagrange Polynoms)
        w1= -( (x1 - x3).*(2*x1 - 3*x2 + x3) ) ./(6*(x1 - x2));
        w2= -(x1 - x3).^3 ./(6*(x1 - x2).*(x2 - x3));
        w3= ((x1 - x3).*(x1 - 3*x2 + 2*x3))./(6*(x2 - x3));
        
        % sums to get single interval integration values
        Inter= w1.*f1 + w2.*f2 + w3.*f3;
        
        % sum over all Intervalls to get final value
        I=sum(Inter);
        

        
        if cum
           % makes shure value of integral is avaivable at every node
           %    -> Standard simpson-law calculates integral over 2 nodes [x_i, x_i+2] 
            
           % Midpoint integral values
           w1M= -((x1 - x2).*(2*x1 + x2 - 3*x3))./(6*(x1 - x3));
           w2M= -((x1 - x2).*(x1 + 2*x2 - 3*x3))./(6*(x2 - x3));
           w3M=  (x1 - x2).^3./(6*(x1 - x3).*(x2 - x3));
           
           InterM = w1M.*f1 + w2M.*f2 + w3M.*f3;
           fint=zeros(size(x));
           fiT=cumsum(Inter);
           fint(3:2:end  )=fiT;
           fint(2:2:end-1)=fint(1:2:end-2) + InterM;
        end
        
        
    else
        % not a even number of intervalls -> uses Midpoint rule for first and simpson for rest
        % First intervall with midpoint rule
        if ~cum
            I= NumInt( f(2:end),x(2:end), 1 ) + 0.5*(x(2)-x(1))*(f(1)+f(2)) ;
        else
           [Itmp,ftmp]=  NumInt( f(2:end),x(2:end), 1, true );
           first=0.5*(x(2)-x(1))*(f(1)+f(2)) ;
           ftmp=ftmp + first;
           I= Itmp + first;
           fint=[0;ftmp];           
        end
    end
    
    
elseif mode==2
    %---------------------------------------------------------
    %           midpoint law
    %---------------------------------------------------------
    
    % mitpoint values of Cf
    fM= 0.5*( f(2:end)+f(1:end-1) );

    % panel length vector
    dx= x(2:end)-x(1:end-1);

    if cum;
        fint= [0;cumsum( fM.*dx )]; 
        I=fint(end);
    else
        I= sum( fM.*dx );
    end
end


end

% %     Herleitung des Verfahrens
% % -------------------------------------------------

% %Koeffizienten für Simpsonregel integrieren
% syms x a b c
% 
% %Lagrange-Funktionen
% L1= (x-b)/(a-b) * (x-c)/(a-c);
% L2= (x-a)/(b-a) * (x-c)/(b-c);
% L3= (x-a)/(c-a) * (x-b)/(c-b);
% 
%  % Analytische Integration er Lagrange-Funktionen über Intervall [x_i, x_i+2] 
% int1= int(L1,x,[a c]);
% int2= int(L2,x,[a c]);
% int3= int(L3,x,[a c]);

%  % Analytische Integration er Lagrange-Funktionen über Intervall [x_i, x_i+1] 
%       -> Für kummulatives integral um Werte an allen Knoten zu haben
% int1M= int(L1,x,[a b]);
% int2M= int(L2,x,[a b]);
% int3M= int(L3,x,[a b]);


% % Beispiel aus Schwarz 2011 zur Verifikation
% %--------------------------------------------
% % analytisch -> Integral ergibt 1
% syms x
% f= 5/(exp(pi)-2)*exp(2*x)*cos(x);
% erg=int(f,x,[0 0.5*pi]);
% 
% ft=@(gg) 5/(exp(pi)-2)*exp(2*gg).*cos(gg);
% xx= 0:pi/8:pi/2; % Schrittweite pi/8
% ff=ft(xx);
% NumInt(ff,xx,1) % mit Mittelpunktsregel
% NumInt(ff,xx,2) % mit Simpson-Regel  


