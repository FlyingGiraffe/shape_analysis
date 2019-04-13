addpath('utils/');

meshfile = 'meshes/166.off';
[X,T] = readOff(meshfile);
nv = size(X,1);
nt = size(T,1);

A = surfaceArea(X,T);
fprintf('The surface area of %s is %f\n', meshfile, A);

% Sanity checks: Laplacian is symmetric and positive definite
L = cotLaplacian(X,T);
[~,p] = chol(L);
fprintf('\nIf %d is 0, then L is PSD\n', p);
fprintf('Symmetry: %d\n', norm(L-L','fro'));

% Divided differences approximation
G = gradApprox(X,T);
fprintf('Difference between gradient and cotan Laplacian %d\n', norm(.5*L*X-G,'fro'));

H = meanCurvature(X,T);

% to determine the value range of colorbar
cmax = max(H);
cmin = min(H);
range = cmax-cmin;
cmax = range/2;
cmin = -range/2;

showDescriptor(X,T,H);
caxis([cmin cmax]);

meshname = strsplit(meshfile,'/');
meshname = strsplit(cell2mat(meshname(end)),'.');
meshname = cell2mat(meshname(1));
saveas(gcf,['4e_mean_curvature_mesh-' meshname '.png']);

%% Mean curvature flow (explicit)
maxiters=1000;
Xt = X;
tau = 0.000001*A;

% ADD CODE FOR EXPLICIT INTEGRATOR HERE %%%%
cotL = cotLaplacian(X,T);
for t=1:maxiters
    Xt = Xt-tau*cotL*Xt./barycentricArea(Xt,T);
    %Xt = Xt-tau*cotLaplacian(Xt,T)*Xt./barycentricArea(Xt,T);
    if mod(t,maxiters/4) == 0
        H = meanCurvature(Xt,T);
        showDescriptor(Xt,T,H);
        caxis([cmin cmax]);
        %saveas(gcf,['5b_mean_curvature_flow_explicit_mesh-' meshname '_step-' int2str(t) '.png']);
        saveas(gcf,['5b_mean_curvature_flow_explicit_nonsing_mesh-' meshname '_step-' int2str(t) '.png']);
    end
end
% END HOMEWORK ASSIGNMENT %%%%
%%  Mean curvature flow (semi-implicit)
maxiters=1000;
Xt = X;
tau = 0.000001*A;

% ADD CODE FOR IMPLICIT INTEGRATOR HERE %%%%
cotL = cotLaplacian(X,T);
for t=1:maxiters
    Xt = (speye(nv)+tau*cotL./barycentricArea(Xt,T))\Xt;
    %Xt = (speye(nv)+tau*cotLaplacian(Xt,T)./barycentricArea(Xt,T))\Xt;
    if mod(t,maxiters/4) == 0
        H = meanCurvature(Xt,T);
        showDescriptor(Xt,T,H);
        caxis([cmin cmax]);
        %saveas(gcf,['5b_mean_curvature_flow_implicit_mesh-' meshname '_step-' int2str(t) '.png']);
        saveas(gcf,['5b_mean_curvature_flow_implicit_nonsing_mesh-' meshname '_step-' int2str(t) '.png']);
    end
end
% END HOMEWORK ASSIGNMENT %%%%

%% Function definitions
% ADD CODE TO COMPUTE SURFACE AREA HERE %%%%%%%%%%
function [A] = surfaceArea(X,T)
    A = sum(arrayfun(@triangleArea,T(:,1),T(:,2),T(:,3)));
    function tarea = triangleArea(T1,T2,T3)
        tarea = norm(cross(X(T2,:)-X(T1,:),X(T3,:)-X(T2,:)))/2;
    end
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%

% ADD CODE TO COMPUTE COTANGENT LAPLACIAN HERE %%%%%%%%%
function [L] = cotLaplacian(X,T)
    nv = size(X,1);
    nt = size(T,1);
    % indices and entries for creating the sparse matrix
    I = zeros(nt*9,1);
    J = zeros(nt*9,1);
    V = zeros(nt*9,1);
    for t = 1:nt
        vertex1 = X(T(t,1),:); vertex2 = X(T(t,2),:); vertex3 = X(T(t,3),:);
        edge12 = vertex2-vertex1;
        edge23 = vertex3-vertex2;
        edge31 = vertex1-vertex3;
        cot1 = -dot(edge12,edge31)/norm(cross(edge12,edge31));
        cot2 = -dot(edge23,edge12)/norm(cross(edge23,edge12));
        cot3 = -dot(edge31,edge23)/norm(cross(edge31,edge23));
        % the diagonal entries
        I(t*9-8) = T(t,1); I(t*9-7) = T(t,2); I(t*9-6) = T(t,3);
        J(t*9-8) = T(t,1); J(t*9-7) = T(t,2); J(t*9-6) = T(t,3);
        V(t*9-8) = (cot2+cot3); V(t*9-7) = (cot1+cot3); V(t*9-6) = (cot1+cot2);
        % the adjacent vertices
        % vertex1 and vertex2 are adjacent
        I(t*9-5) = T(t,1); I(t*9-4) = T(t,2);
        J(t*9-5) = T(t,2); J(t*9-4) = T(t,1);
        V(t*9-5) = -cot3; V(t*9-4) = -cot3;
        % vertex2 and vertex3 are adjacent
        I(t*9-3) = T(t,2); I(t*9-2) = T(t,3);
        J(t*9-3) = T(t,3); J(t*9-2) = T(t,2);
        V(t*9-3) = -cot1; V(t*9-2) = -cot1;
        % vertex1 and vertex3 are adjacent
        I(t*9-1) = T(t,1); I(t*9) = T(t,3);
        J(t*9-1) = T(t,3); J(t*9) = T(t,1);
        V(t*9-1) = -cot2; V(t*9) = -cot2;
    end
    % recurring entries will be summed up
    L = sparse(I,J,V,nv,nv);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%

% ADD CODE TO COMPUTE DIVIDED DIFFERENCES APPROXIMATION HERE %%%%%
function [G] = gradApprox(X,T)
    nv = size(X,1);
    nt = size(T,1);
    G = zeros(nv,3);
    h = 0.001;
    for axis = 1:3
        % the change of area w.r.t each point is local
        for t = 1:nt
            v1 = X(T(t,1),:); v2 = X(T(t,2),:); v3 = X(T(t,3),:);
            e = zeros(1,3);
            e(axis) = 1;
            G(T(t,1),axis) = G(T(t,1),axis)+(tarea(v1+e*h,v2,v3)-tarea(v1-e*h,v2,v3))/(2*h);
            G(T(t,2),axis) = G(T(t,2),axis)+(tarea(v1,v2+e*h,v3)-tarea(v1,v2-e*h,v3))/(2*h);
            G(T(t,3),axis) = G(T(t,3),axis)+(tarea(v1,v2,v3+e*h)-tarea(v1,v2,v3-e*h))/(2*h);
        end
    end
    function tarea = tarea(X1,X2,X3)
        tarea = norm(cross(X1-X2,X3-X2))/2;
    end
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%%%%

% ADD CODE TO COMPUTE THE BARYCENTRIC AREA VECTOR HERE %%%%%%%%%
function [M] = barycentricArea(X,T)
    nv = size(X,1);
    nt = size(T,1);
    M = zeros(nv,1);
    for t = 1:nt
        % the barycentric area of each vertex in triangle t
        A = norm(cross(X(T(t,2),:)-X(T(t,1),:),X(T(t,3),:)-X(T(t,2),:)))/2;
        M(T(t,1)) = M(T(t,1))+A/3;
        M(T(t,2)) = M(T(t,2))+A/3;
        M(T(t,3)) = M(T(t,3))+A/3;
    end
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%%%

% ADD CODE TO COMPUTE POINTWISE MEAN CURVATURE HERE %%%%%%%%%%%
function [H] = meanCurvature(X,T)
    Hn = 0.5*cotLaplacian(X,T)*X;
    H = sqrt(sum(Hn.*Hn,2))./barycentricArea(X,T);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%