function blo=addBlowing(blo,prf)
%ADDBLOWING determines the blowing velocity for each nodes


blo.Vb=zeros(prf.N,1);
% if blowing inactive -> return
if ~blo.active; return; end


x = {prf.nodes.X(1:prf.M), prf.nodes.X(prf.M:end)};


% if ~isfield(blo,'ArcLengthMode')
%     blo.ArcLengthMode=false;
% end



if blo.Mode==1 || blo.Mode==2  
    
    % find blowing spots and set corresponding velocities
    
    for side=1:2 % loop of sides
        for i=1:numel(blo.x{side})
            if blo.Mode==1 
                % x_start and L for both sides
                ind=find( x{side} > blo.x{side}(i) ...
                        & x{side}<(blo.x{side}(i)+blo.L{side}(i)) );
            elseif blo.Mode==2
                % x_mid and L for both sides
                ind = find( x{side}>(blo.x{side}(i)-blo.L{side}(i)/2)...
                          & x{side}<(blo.x{side}(i)+blo.L{side}(i)/2)); 
            end 
                 blo.Vb(ind + (side==2)*(prf.M-1)) = blo.A{side}(i);
        end
    end
    
elseif blo.Mode==3
% mode for arc length
      
    % error for wrong parameters
    if length(blo.s_change)~=length(blo.NewA)
        error('The Vectors blo.NewA and blo.s_change musst have the same length -> check the parameters in the script')
    elseif ~isempty(find(blo.s_change>=1,1))
        error('The elements of blo.s_change musst be between 0 and 1 -> check the parameters in the script')
    end
    s_ch=blo.s_change*prf.s(end);    


    %ind=[]; % vector with corresponding indices 
    for i=1:length(s_ch)
        %ind=[ind, find( prf.s>s_ch(i),1,'first' ) ];
        %blo.Vb(ind(i):prf.N)=blo.NewA(i);

        % find nodes
        ind=find( prf.s>s_ch(i),1,'first' );
        % change blowing / suction velocity
        blo.Vb(ind:prf.N)=blo.NewA(i);
    end


end

end

