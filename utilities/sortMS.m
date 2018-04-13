function out = sortMS( dat)
%SORTMS sorts for upper and lower airfoil part
%       dat= [x,y,.....]
%       -> out: 1. upper part from x=1,..,0
%               2. lower part from x=0,..,1

%split lower and upper part of airfoil
indU=find( dat(:,2)>0  ); % upper
indL=find( dat(:,2)<0  ); % lower

tmpU=dat(indU,:);
tmpL=dat(indL,:);

% sorts after x
[~,sortU]=sort(tmpU(:,1),'descend');
tmpU=tmpU(sortU,:);

[~,sortL]=sort(tmpL(:,1));
tmpL=tmpL(sortL,:);

out= [tmpU;tmpL];



end

