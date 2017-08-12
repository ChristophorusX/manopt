function [clustering, eigvalues, V] =  spectral_clustering_rw(W, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function only dealing with L_{rw}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,~]=size(W);
V = zeros(n,k,2);

D2 = diag(1./sum(W,2));
L2 = speye(n) - D2*W;
[X2, ~] = svd(L2);
X2 = X2(:,end-k+1:end);
rng(2);
[grp2,~] = kmeans(X2, k);
clustering = 2*(grp2 - 1) - 1;
eigvalues(:,2)=eig(L2);
V(:,:,2)=X2;
%

end
