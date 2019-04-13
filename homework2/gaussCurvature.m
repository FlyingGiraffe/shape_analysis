addpath('utils/');

meshfile = 'meshes/moomoo.off';
[X,T] = readOff(meshfile);
nv = size(X,1);
nt = size(T,1);

% FILL IN TWO METHODS FOR COMPUTING GAUSSIAN CURVATURE HERE %%%%%
% COMPUTE GAUSSIAN USING SECOND FUNDAMENTAL FORM
% compute the normal vector at each vertex
vnorm = zeros(nv,3); % normal of each vertex
for t = 1:nt
    vertex1 = X(T(t,1),:); vertex2 = X(T(t,2),:); vertex3 = X(T(t,3),:);
    edge12 = vertex2-vertex1;
    edge23 = vertex3-vertex2;
    edge31 = vertex1-vertex3;
    tnorm = cross(edge12,edge23)/norm(cross(edge12,edge23));
    A = norm(cross(edge12,edge23))/2;
    vnorm(T(t,1),:) = vnorm(T(t,1),:)+tnorm*A;
    vnorm(T(t,2),:) = vnorm(T(t,2),:)+tnorm*A;
    vnorm(T(t,3),:) = vnorm(T(t,3),:)+tnorm*A;
end
vnorm = vnorm./sqrt(sum(vnorm.*vnorm,2)); % normalize
% compute the second fundamental form at each vertex
Mv = zeros(nv*3,3);
v_total_weight = zeros(nv,1);
for t = 1:nt
    vertex1 = X(T(t,1),:); vertex2 = X(T(t,2),:); vertex3 = X(T(t,3),:);
    edge12 = vertex2-vertex1;
    edge23 = vertex3-vertex2;
    edge31 = vertex1-vertex3;
    % vertex 1
    t12 = (eye(3)-vnorm(T(t,1),:).'*vnorm(T(t,1),:))*edge12.'; t12 = t12/norm(t12);
    t13 = (eye(3)-vnorm(T(t,1),:).'*vnorm(T(t,1),:))*(-edge31.'); t13 = t13/norm(t13);
    k12 = 2*dot(vnorm(T(t,1),:),edge12)/sum(edge12.*edge12);
    k13 = 2*dot(vnorm(T(t,1),:),-edge31)/sum(edge31.*edge31);
    % vertex 2
    t21 = (eye(3)-vnorm(T(t,2),:).'*vnorm(T(t,2),:))*(-edge12.'); t21 = t21/norm(t21);
    t23 = (eye(3)-vnorm(T(t,2),:).'*vnorm(T(t,2),:))*edge23.'; t23 = t23/norm(t23);
    k21 = 2*dot(vnorm(T(t,2),:),-edge12)/sum(edge12.*edge12);
    k23 = 2*dot(vnorm(T(t,2),:),edge23)/sum(edge23.*edge23);
    % vertex 3
    t31 = (eye(3)-vnorm(T(t,3),:).'*vnorm(T(t,3),:))*edge31.'; t31 = t31/norm(t31);
    t32 = (eye(3)-vnorm(T(t,3),:).'*vnorm(T(t,3),:))*(-edge23.'); t32 = t32/norm(t32);
    k31 = 2*dot(vnorm(T(t,3),:),edge31)/sum(edge31.*edge31);
    k32 = 2*dot(vnorm(T(t,3),:),-edge23)/sum(edge23.*edge23);
    % area of triangle t = weight on edge uv in this triangle
    A = norm(cross(edge12,edge23))/2;
    Mv(T(t,1)*3-2:T(t,1)*3,:) = Mv(T(t,1)*3-2:T(t,1)*3,:) + A*k12*(t12*t12.') + A*k13*(t13*t13.');
    Mv(T(t,2)*3-2:T(t,2)*3,:) = Mv(T(t,2)*3-2:T(t,2)*3,:) + A*k21*(t21*t21.') + A*k23*(t23*t23.');
    Mv(T(t,3)*3-2:T(t,3)*3,:) = Mv(T(t,3)*3-2:T(t,3)*3,:) + A*k31*(t31*t31.') + A*k32*(t32*t32.');
    % the weight w should be normalized to sum to 1 at last
    v_total_weight(T(t,1)) = v_total_weight(T(t,1))+A*2;
    v_total_weight(T(t,2)) = v_total_weight(T(t,2))+A*2;
    v_total_weight(T(t,3)) = v_total_weight(T(t,3))+A*2;
end
gaussianCurvature1 = zeros(nv,1);
for v = 1:nv
    II = Mv(v*3-2:v*3,:)/v_total_weight(v);
    lambda = eigs(II,2); % the two non-zero eigenvalues
    gaussianCurvature1(v) = (3*lambda(1)-lambda(2))*(3*lambda(2)-lambda(1));
end

% STRUCTURE-PRESERVING GAUSSIAN CURVATURE
angle = zeros(nv,1);
varea = zeros(nv,1);
for t = 1:nt
    vertex1 = X(T(t,1),:); vertex2 = X(T(t,2),:); vertex3 = X(T(t,3),:);
    edge12 = vertex2-vertex1;
    edge23 = vertex3-vertex2;
    edge31 = vertex1-vertex3;
    % the interior angle at each vertex in triangle t
    theta1 = acos(dot(edge12,-edge31)/(norm(edge12)*norm(edge31)));
    theta2 = acos(dot(edge23,-edge12)/(norm(edge23)*norm(edge12)));
    theta3 = acos(dot(edge31,-edge23)/(norm(edge31)*norm(edge23)));
    angle(T(t,1)) = angle(T(t,1))+theta1;
    angle(T(t,2)) = angle(T(t,2))+theta2;
    angle(T(t,3)) = angle(T(t,3))+theta3;
    % the area of Voronoi cell at each vertex in triangle t
    A = norm(cross(edge12,edge23))/2;
    varea(T(t,1)) = varea(T(t,1))+A/3;
    varea(T(t,2)) = varea(T(t,2))+A/3;
    varea(T(t,3)) = varea(T(t,3))+A/3;
end
gaussianCurvature2 = (2*pi-angle)./varea;
% END HOMEWORK PROBLEM %%%%%%%%%%%%%%%%%%%%

meshname = strsplit(meshfile,'/');
meshname = strsplit(cell2mat(meshname(end)),'.');
meshname = cell2mat(meshname(1));

% to determine the value range of colorbar
cmax = max([gaussianCurvature1;gaussianCurvature2]);
cmin = min([gaussianCurvature1;gaussianCurvature2]);
range = cmax-cmin;
cmax = range/7;
cmin = -range/7;

showDescriptor(X, T, gaussianCurvature1);
caxis([cmin cmax]);
saveas(gcf,['3_gauss_curvature1_mesh-' meshname '.png']);
showDescriptor(X, T, gaussianCurvature2);
caxis([cmin cmax]);
saveas(gcf,['3_gauss_curvature2_mesh-' meshname '.png']);