function blo=addBlowing(blo,prf)



if ~isfield(blo,'ArcLengthMode')
    blo.ArcLengthMode=false;
end

if blo.active
    x = {prf.nodes.X(1:prf.M), prf.nodes.X(prf.M:end)};
    if blo.ArcLengthMode==false
        %------------------ first mode ---------------------
        for side=1:2

            for i=1:numel(blo.x{side})

                ind = find( x{side}>(blo.x{side}(i)-blo.L{side}(i)/2)...
                          & x{side}<(blo.x{side}(i)+blo.L{side}(i)/2));

                blo.Vb(ind + (side==2)*(prf.M-1)) = blo.A{side}(i);

            end

        end
        
    else
        %------------------ second mode --------------------
        
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

