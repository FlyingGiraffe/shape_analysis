addpath('utils/');

filename = 'swissroll.txt';
data = readtable(filename);

X = [data.Var1 data.Var2 data.Var3];
n = size(X,1);
S = 10*ones(n,1);
C = X(:,3);
figure;
scatter3(X(:,1),X(:,2),X(:,3),S,C,'filled');
axis equal;

%%% ADD CODE FOR YOUR TWO METHODS HERE %%%

%% SMACOF
X_SMACOF = rand(2,n);
D0 = zeros(n);
for i=1:n
    for j=1:n
        D0(i,j) = norm(X(i,:)-X(j,:));
    end
end
maxStep = 100;
for k = 1:maxStep
    X_SMACOF = X_SMACOF * B(X_SMACOF,D0) * (eye(n)-ones(n)./n) ./(2*n);
end
figure;
scatter(X_SMACOF(1,:),X_SMACOF(2,:),S,C,'filled')
axis equal;

%% K-ISOMAP
D = distance(X,10);
X_K_ISO = MDS(D);

figure;
scatter(X_K-ISO(1,:),X_K-ISO(2,:),S,C,'filled')
axis equal;

%% eps-ISOMAP
eps = 6;
D_eps = distance_eps(D0,eps);
X_eps_ISO = MDS(D_eps);

figure;
scatter(X_eps_ISO(1,:),X_eps_ISO(2,:),S,C,'filled')
axis equal;

%% Function definitions

%%% SMACOF %%%
% Compute matrix B(X) at each iteration step
function [matB] = B(X,D0)
    n = size(X,2);
    matB = zeros(n);
    for i=1:n
        for j=1:n
            if i~=j && norm(X(:,i)-X(:,j))~=0
                matB(i,j) = -D0(i,j)/norm(X(:,i)-X(:,j));
            end  
        end
    end
    for i=1:n
        matB(i,i) = -sum(matB(i,:));
    end
end
%%% END SMACOF %%%

%%% K-ISOMAP %%%
% Find the K nearest neighbours of each point, compute the distance matrix using Floyd-Warshall algorithm
function [D] = distance(X,K)
    n = size(X,1);
    D = Inf(n);
    [Idx,D_edge] = knnsearch(X,X,'K',K+1);
    for i=1:n
        D(i,i) = 0;
    end
    for i=1:n
        for k=1:K
            D(i,Idx(i,k)) = D_edge(i,k);
            D(Idx(i,k),i) = D_edge(i,k);
        end
    end
    for k=1:n
        for i=1:n
            for j=1:n
                if D(i,j) > D(i,k)+D(k,j)
                    D(i,j) = D(i,k)+D(k,j);
                end
            end
        end
    end
end

function [D] = distance_eps(D0,eps)
    n = size(D0,1);
    D = D0;
    D(D>eps) = inf;
    for k=1:n
        for i=1:n
            for j=1:n
                if D(i,j) > D(i,k)+D(k,j)
                    D(i,j) = D(i,k)+D(k,j);
                end
            end
        end
    end
end

% Classical MDS
function [X_MDS] = MDS(D)
    n = size(D,1);
    P = D.*D;
    G = -(eye(n)-ones(n)./n)*P*(eye(n)-ones(n)./n)./2;
    [V,Diag] = eigs(G,2);
    X_MDS = sqrt(Diag)*V.';
end
%%% END K-ISOMAP %%%

%%% END HOMEWORK ASSIGNMENT %%%