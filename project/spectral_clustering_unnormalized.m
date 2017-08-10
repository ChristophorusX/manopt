function [grp, eigvalues, V] =  spectral_clustering_unnormalized(W, k)

[n,~]=size(W);
V = zeros(n,k,3);

D1 = diag(sum(W,2));
[X1, ~]=svd(D1-W);
X1 = X1(:,end-k+1:end);
rng(2);
[grp1,~] = kmeans(X1,k);
eigvalues(:,1)=eig(D1-W);
V(:,:,1)=X1;
%
D2 = diag(1./sum(W,2));
L2 = speye(n) - D2*W;
[X2, ~] = svd(L2);
X2 = X2(:,end-k+1:end);
rng(2);
[grp2,~] = kmeans(X2, k);
eigvalues(:,2)=eig(L2);
V(:,:,2)=X2;
%
D3 = diag(1 ./ sqrt(sum(W, 2)));
L3 = speye(n)- D3 * W * D3;
[X3, ~] = svd(L3);
X3 = X3(:,end-k+1:end);
Y = X3 ./ repmat(sqrt(sum(X3.^2, 2)), 1, k);
rng(2);
% [grp3,~] = kmeans(Y, k);
grp3 = zeros(n,1);
eigvalues(:,3)=eig(L3);
V(:,:,3)=Y;
grp = grp1;

end
