function [A] = massMatrix(X,T)

nv = size(X,1);
% Triangle areas
N = cross(X(T(:,1),:)-X(T(:,2),:), X(T(:,1),:)-X(T(:,3),:));
Ar = .5*normv(N);

% Vertex areas = sum triangles nearby
I = [T(:,1);T(:,2);T(:,3)];
S = [Ar;Ar;Ar]/3;
A = sparse(I,I,S,nv,nv);

function nn = normv(V)
nn = sqrt(sum(V.^2,2));

