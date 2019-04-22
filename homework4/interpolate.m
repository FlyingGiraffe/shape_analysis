clear;
addpath('utils/');
%% Load mesh
filename1 = 'meshes/mushroom.off';
[X1,T1] = readOff(filename1);
X1 = rescale(X1, .03, .97);
mesh1 = getMeshData(X1,T1);

filename2 = 'meshes/duck.off';
[X2,T2] = readOff(filename2);
X2 = rescale(X2, .03, .97);
mesh2 = getMeshData(X2,T2);

%% Convert the meshes into voxels
N = 100;
t = linspace(0,1,N);
volume1 = inpolyhedron(T1,X1,t,t,t);
volume2 = inpolyhedron(T2,X2,t,t,t);
%% View the voxels
opts.alpha = 1; % transparency
opts.color = [1 0 0];
figure;
plot_isosurface(volume1,opts);

opts.color = [0 0 1];
figure;
plot_isosurface(volume2,opts);

%% Interpolation
volume1 = volume1./sum(volume1,'all');
volume2 = volume2./sum(volume2,'all');
time = 1;
kernel = @(x) imgaussian(x,time,N);
H0 = max(entropy(volume1),entropy(volume2));
niter = 1000;
n_interp = 3;
interp = {};
for i = 1:n_interp
    alpha1 = i/(n_interp+1); alpha2 = 1-alpha1;
    interp{i} = barycenter(volume1,volume2,alpha1,alpha2,kernel,H0,niter);
end

%% Convert the interpolated volumes into meshes
inter_meshes = {}
for i = 1:n_interp
    figure
    inter_meshes{i} = plot_isosurface(interp{i},opts);
end

%% Plot the results
figure
opts.color = [0 1 0];
plot_mesh(X1,T1,opts);
for i = 1:n_interp
    figure
    alpha1 = i/(n_interp+1); alpha2 = 1-alpha1;
    opts.color = [0 alpha1 alpha2];
    plot_mesh(inter_meshes{i}.vertices,inter_meshes{i}.faces,opts);
end
figure
opts.color = [0 0 1];
plot_mesh(X2,T2,opts);