function [groups, eigvalues, V] =  SpectralClustering(W, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to implement the 3 variations of the Specral Clustering
% Input
% W: n * n adjacency matrix
% k: the number of desired clusters
% Output
% groups : n*3 vector, indicating the cluster assigment
% === first column indicates the cluster assigment based on SC1(unnormalized)
% === second column indicates the cluster assigment based on SC2(rw-noramlized)
% === third column indicates the cluster assigment based on SC3(sym-normalized)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,~] = size(W);
V = zeros(n,k,3);

D_1 = diag(sum(W,2));
[X_1, ~] = svd(D_1-W);
X_1 = X_1(:,end-k+1:end);
rng(2);
[group_1,~] = kmeans(X_1,k);
eigvalues(:,1) = eig(D_1-W);
V(:,:,1) = X_1;
%
D_2 = diag(1./sum(W,2));
L_2 = speye(n) - D_2*W;
[X_2, ~] = svd(L_2);
X_2 = X_2(:,end-k+1:end);
rng(2);
[group_2,~] = kmeans(X_2, k);
eigvalues(:,2) = eig(L_2);
V(:,:,2) = X_2;
%
D_3 = diag(1 ./ sqrt(sum(W, 2)));
L_3 = speye(n)- D_3 * W * D_3;
[X_3, ~] = svd(L_3);
X_3 = X_3(:,end-k+1:end);
Y = X_3 ./ repmat(sqrt(sum(X_3.^2, 2)), 1, k);
rng(2);
[group_3,~] = kmeans(Y, k);
eigvalues(:,3) = eig(L_3);
V(:,:,3) = Y;
groups = [group_1,group_2,group_3];
