clear;
addpath('utils/');
%% Load mesh
filename = 'meshes/moomoo.off';
[X,T] = readOff(filename);
nv = size(X,1);
nt = size(T,1);
mesh = getMeshData(X,T);

%% Define distributions
p = zeros(nv,1);
source = 5;
p(source) = 1;

%q = 1/nv*ones(nv,1);
%q = zeros(nv,1);
%target = 1;
%q(target) = 1;

showDescriptor(mesh,p);
%showDescriptor(mesh,q);

%% Regularized EMD code
time = .00001;
%time = 0.05;
steps = 5;
%kernel = @(x) heatDiffusion(x,mesh,time,steps);

%alpha = 0.1;
alpha = .00002;
niter = 1000;

% YOUR CODE HERE TO PERFORM SINKHORN ITERATIONS %%%
W = zeros(nv,1);
for target = 1:nv
    q = zeros(nv,1);
    q(target) = 1;
    W(target) = EMD(p,q,mesh,time,steps,alpha,niter);
end
showDescriptor(mesh,W);
% END ASSIGNMENT %%%

%% Function definitions
function W = EMD(p,q,mesh,time,steps,alpha,niter)
    kernel = @(x) heatDiffusion(x,mesh,time,steps);
    nv = size(mesh.vertices,1);
    v = ones(nv,1); w = ones(nv,1);
    threshold = 1e-6;
    W_old = 0;
    for k = 1:niter
        v = p./kernel(w);
        w = q./kernel(v);
        logv = log(v); logv(v<=0) = 0;
        logw = log(w); logw(w<=0) = 0;
        W = alpha*sum(p.*logv+q.*logw);
        if norm(W_old-W) < threshold  &&  k > 2
            break;
        end
        W_old = W;
    end
    % If an entry of v(resp. w) is 0, the corresponding entry in p(resp. q) is also 0
    logv = log(v); logv(v<=0) = 0;
    logw = log(w); logw(w<=0) = 0;
    W = alpha*sum(p.*logv+q.*logw);
end