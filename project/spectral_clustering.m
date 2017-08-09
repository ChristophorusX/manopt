function [ clustering, eigenvalues, eigenvectors ] = spectral_clustering( G )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform spectral clustering function on the graph, trying to recover the
% planted partition by grouping the eigenvectors of the adjacency matrix
% with k-means. (The helping function is in seperate function file.)
% It receives
% - A matrix G to perform spectral clustering on.
% It returns
% - The resulting clustering of all the vertices in 0 and 1.
% - All the eigenvalues of the matrix.
% - All the eigenvectors of the matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[groups,eigenvalues,eigenvectors] = SpectralClustering(G,2); % Seperate into 2 groups
clustering = groups(:,3)-1;
% display(clustering)

end

