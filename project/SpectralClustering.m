function [grps, eigvalues, V] =  SpectralClustering(W, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to implement the 3 variations of the Specral Clustering
% Input 
% W: n * n adjacency matrix
% k: the number of desired clusters
% Output
% grps : n*3 vector, indicating the cluster assigment 
% === first column indicates the cluster assigment based on SC1(unnormalized) 
% === second column indicates the cluster assigment based on SC2(rw-noramlized) 
% === third column indicates the cluster assigment based on SC3(sym-normalized) 
%
% EC 500: Learning from Data
% Fall Semester, 2015
% Instructor: Prakash Ishwar
% MATLAB Exercise 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(W);
V = zeros(n,k,3);

D1 = diag(sum(W,2));
[X1, d1]=svd(D1-W);
X1 = X1(:,end-k+1:end);
rng(2);
[grp1,~] = kmeans(X1,k);
eigvalues(:,1)=eig(D1-W);
V(:,:,1)=X1;
% 
D2 = diag(1./sum(W,2));
L2 = speye(n) - D2*W;
[X2, d2] = svd(L2);
X2 = X2(:,end-k+1:end);
rng(2);
[grp2,~] = kmeans(X2, k);
eigvalues(:,2)=eig(L2);
V(:,:,2)=X2;
%
D3 = diag(1 ./ sqrt(sum(W, 2)));
L3 = speye(n)- D3 * W * D3;
[X3, d3] = svd(L3);
X3 = X3(:,end-k+1:end);
Y = X3 ./ repmat(sqrt(sum(X3.^2, 2)), 1, k);
rng(2);
[grp3,~] = kmeans(Y, k);
eigvalues(:,3)=eig(L3);
V(:,:,3)=Y;
grps = [grp1,grp2,grp3];
