function x = BackwardsSubstitution( A, b,ind)
%BACKWARDSSUBSTITUTION führt eine Rücksubstitution durch. A ist in LR-Form,
%                       b die rechte Seite. Ind ist ein Vektor, der die
%                       Zeilenvertauschungen enthält

n=length(b);
bper=b(ind);
c=bper;

% rechte Seite für LR-Form
for i=2:n
    c(i)=c(i) - sum(A(i,1:i-1)*c(1:i-1));
end

% Rücksubstitutuin
%x=zeros(n,1);
x=zeros(size(b));

for i=n:-1:1
   s=c(i);
   if i<n % new rhs
    s=s-sum(A(i,i+1:n)*x(i+1:n));
   end
   x(i)=s/A(i,i);
end


end

