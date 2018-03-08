function out = Grading( in, xEND, N ,getGrad)
%GRADING if getGrad: gets the grading for a fixed number of point, a fixed starting length and a fixed end length
%               -> in is starting step size
%        if not  getGrad: gets starting step size with for a fixed number of point, a fixed grading and a fixed end length
      
if nargin==3
    getGrad=true;
end

if getGrad
    dx1=in;
    
    Nq=xEND/dx1;
    
    % quadratic initial guess
    a= (N-1)*(N-2)*(N-3)/6;
    b= (N-1)*(N-2)/2;
    c= (N-1)-Nq;
    
    d= max(0,b^2-4*a*c);
    
    arg=1/(N-1);
    if N==3
       grad= -c/b+1;
    else
       grad=(-b+sqrt(d)) /(2*a) + 1;
    end
    
    k=0; dgr=1;
    while abs(dgr)>1e-5 && k<60
        Nnew= (grad^(N-1)-1)/(grad-1);
        res= Nnew^arg- Nq^arg ;
        
        dres_dR= arg*Nnew^arg*( (N-1)*grad^(N-1)-Nnew)/(grad^(N-1)-1);
        dgr=-res/dres_dR;
        grad=grad + dgr;
        
        k=k+1;
    end
    
    out=grad;

elseif~getGrad 
    grad=in;
    ind=1:N-2;
    gr=grad.^ind;
    dx1= xEND/(1+sum( gr ));
    out=dx1;
end

% ind=1:N-2;
% gr=grad.^ind;
% x= cumsum( dx1*[0,1,gr] )';

end

