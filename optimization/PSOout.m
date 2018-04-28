function stop = PSOout( optimValues,state )
%PSOOUT

stop=false;
if optimValues.iteration>38
   stop=true; 
end
x= optimValues.swarm;
y=-optimValues.swarmfvals;
global PSOerg
PSOerg= [PSOerg;x,y];

end

