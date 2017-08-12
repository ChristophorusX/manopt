function [clustering, eigvalues, V] =  spectral_clustering_rw(W, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function only dealing with L_{rw}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,~] = size(W);
V = zeros(n,k);

D_2 = diag(1./sum(W,2));
L_2 = speye(n) - D_2*W;
[X_2, ~] = svd(L_2);
X_2 = X_2(:,end-k+1:end);
rng(2);
[group_2,~] = kmeans(X_2, k);
clustering = 2*(group_2 - 1) - 1;
eigvalues = eig(L_2);
V(:,:) = X_2;

end
