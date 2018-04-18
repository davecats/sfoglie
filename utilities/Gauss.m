function erg = Gauss( A, b)
%GAUSS Gauß-Algorithmus mit relativer Spaltenmaximumstrategie 
%       -> betragsmäßiges Maximum einer Spalte als Pivotelement
%      gibt rechte und linke Dreiecksmatrix zurück
%      Falls rechte Seite Übergeben -> löst LGS, sonst nur LR-Zerlegung 

n=length(A(:,1));
dt=1;
ind=(1:n)';

for k=1:n-1 % Schleife über Zeilen
   mx=0;
   pk=0;
   
   % Pivot-Element finden und ggf Zeilen tauschen
   for i=k:n 
       s=sum( abs(A(i,k:n)) );
       q=abs(A(i,k))/s;
       if q>mx
         mx=q;
         pk=i;              
       end
   end
   
   if mx==0; disp('gg'); break; end
   if pk~=k
      dt=-dt;
      % Zeile tauschen
      h=A(k,:);
      A(k,:)= A(pk,:);
      A(pk,:)=h;  
      tmp=ind(k);
      ind(k)=ind(pk);
      ind(pk)=tmp;
   end
   %------------------------
   
   dt=dt*A(k,k);
   % Spalte durch aktuelles kk Element teilen
   A(k+1:n,k)=A(k+1:n,k)/A(k,k);
   
   % aktuelle Pivotzeile abziehen
   
   A(k+1:n, k+1:n)= A(k+1:n, k+1:n) - A(k+1:n, k)*A(k,k+1:n);
  
end

erg.dt=dt*A(n,n); % Determinante
erg.L=tril(A,-1); % -1 ohne Hauptdiagonale
erg.L=erg.L+diag(ones(n,1));
erg.R=triu(A);
erg.A=A;
erg.ind=ind;

if nargin==2
    % permutierte rechte Seite
    bper=b(ind);
    c=bper;

    % rechte Seite
    for i=2:n
        c(i)=c(i) - sum(A(i,1:i-1)*c(1:i-1));
    end

    % Rücksubstitutuin
    x=zeros(n,1);

    for i=n:-1:1
       s=c(i);
       if i<n % new rhs
        s=s-sum(A(i,i+1:n)*x(i+1:n));
       end
       x(i)=s/A(i,i);
    end


    erg.c=c;
    erg.x=x;
end



end

