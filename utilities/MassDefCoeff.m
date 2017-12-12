function [D] = MassDefCoeff(B, s,verf)
%MASSDEFCOEFF   calculates the coefficients of the mass defect m_i from B_ij q_j            
%               substituting q_j=(m_j+1-m_j)/(s_j+1-s_j) gives D_ij m_j 
%               default: using forward differences
%               if verf is committed, using central differences

n=length(B(:,1));

if nargin==2
    %forward difference; backward difference for the end node
    tmp=1./(s(2:end)-s(1:end-1)); 
    tmp=ones(n,1)*tmp;
    
    C1=[tmp(:,1:end).*B(:,1:end-1), zeros(n,1)];
    C2=[zeros(n,1), tmp.*B(:,1:end-1)];
    D=C2-C1;   % all terms of forward differencing q_j=2,..,n-1
               %    -> D_ij=B_ij-1/(s_j-s_j-1) - B_ij/(s_j+1-s_j) , j=2...n-1
               
    % last source contribution with backwards diff 
    % -> XFoil: q_N+Nw = 0 -> weglassen
   % last= tmp(:,end).*B(:,end);
   % D(:,end-1:end)=D(:,end-1:end)+[-last, last];

    
    %last= -tmp(:,end).*B(:,end-1);
    %D(:,end-1)=D(:,end-1)+last;
else
    %central difference; using forward diff for the start and backward diff for the end node
    tmp=1./(s(3:end)-s(1:end-2));
    tmp=ones(n,1)*tmp;
    
    C1=[tmp.*B(:,2:end-1),zeros(n,1),zeros(n,1)];
    C2=[zeros(n,1),zeros(n,1),tmp.*B(:,2:end-1)];
    D=C2-C1; % all terms of central differences q_j=2,...,n-1
             %              -> q_j=m_j+1-m_j-1/(s_j+1-s_j-1)
             %              -> D_ij= B_ij-1/(s_j-s_j-2)-B_ij+1/(s_j+2-s_j)
    
    %corrections for start and end point
    first=1/(s(2)-s(1))*B(:,1);
    D(:,1:2)=D(:,1:2) + [-first, first]; % forward diff
    %last=1/(s(end)-s(end-1))*B(:,end);
    %D(:,end-1:end)=D(:,end-1:end) + [-last, last]; % backward diff
end


end

