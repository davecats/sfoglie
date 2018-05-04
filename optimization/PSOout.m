function stop = PSOout( optimValues,state )
%PSOOUT outputfunction for particle swarm optimization 

maxIter=38; % maximum number of iterations
stop=false;
if optimValues.iteration>maxIter
   stop=true; 
end
x= optimValues.swarm;
y=-optimValues.swarmfvals;
global PSOerg
PSOerg= [PSOerg;x,y];

end

