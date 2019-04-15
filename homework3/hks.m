filename = 'human_coarse.off';
[X,T] = readOff(filename);
nv = size(X,1);
nt = size(T,1);
mesh = getMeshData(X,T);

%% Compute 10 eigenfunctions and display eigenfunction 4
k = 10;
[vals, vecs] = laplacianSpectrum(X,T,k);
showDescriptor(mesh,vecs(:,4));

%% Plot HKS on hands and knees of mesh
k = 300;
nSamples = 500;

[vals,vecs] = laplacianSpectrum(X,T,k);

% scale the eigenvectors s.t. their square integral on the surface =1
A = massMatrix(X,T);
square_integral = (abs(vecs).^2).' * diag(A);
vecs = vecs./(square_integral.');

heatSignature = HKS(vals,vecs,nSamples);

% scale by area
heatSignature = heatSignature./sum(A*heatSignature);

% hand1 = 259;
% hand2 = 135;
% knee1 = 232;
% knee2 = 257;

hand1 = 410;
hand2 = 429;
knee1 = 225;
knee2 = 226;

x = 1:nSamples;
f = figure;
subplot(2,2,1);
scatter(x,heatSignature(hand1,:),5);
title('hand 1')
subplot(2,2,2);
scatter(x,heatSignature(hand2,:),5);
title('hand 2')
subplot(2,2,3);
scatter(x,heatSignature(knee1,:),5);
title('knee 1')
subplot(2,2,4);
scatter(x,heatSignature(knee2,:),5);
title('knee 2')

%% COMPUTE DIFFERENCE FUNCTION |HKS(x0) - HKS(x)|_2 HERE ###
x0 = hand2;
result = sum(abs(heatSignature(x0,:)-heatSignature),2);
showDescriptor(mesh,result);

%% COMPLETE THE FUNCTIONS BELOW
function [eigenvalues,eigenvectors] = laplacianSpectrum(X,T,k) 
[W,A] = cotLaplacian(X,T);
[eigenvectors,D] = eigs(W./A,k,'smallestabs');
eigenvalues = diag(D);
end

function heatSignature = HKS(eigenvalues,eigenvectors,nSamples)
% neig = size(eigenvalues,1);

tmin = 4*log(10)/eigenvalues(end);
tmax = 4*log(10)/eigenvalues(2);
ts = logspace(log10(tmin),log10(tmax),nSamples);

% heatSignature = zeros(size(eigenvectors,1),nSamples);
heatSignature = abs(eigenvectors).^2 * exp((-abs(eigenvalues))*ts);
end