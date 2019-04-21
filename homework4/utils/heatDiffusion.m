function result = heatDiffusion(signal,M,time,steps)

h = time/steps;
nv = M.numVertices;

blurInverse = spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;

result = signal;
for i=1:steps
    result = blurInverse \ bsxfun(@times,result,M.areaWeights);
end
