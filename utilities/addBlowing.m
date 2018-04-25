function blo=addBlowing(blo,prf)

x = {prf.nodes.X(1:prf.M), prf.nodes.X(prf.M:end)};

if blo.active
for side=1:2
    
    for i=1:numel(blo.x{side})
        
        ind = find( x{side}>(blo.x{side}(i)-blo.L{side}(i)/2)...
                  & x{side}<(blo.x{side}(i)+blo.L{side}(i)/2));
        
        blo.Vb(ind + (side==2)*(prf.M-1)) = blo.A{side}(i);
        
    end
    
end
end

