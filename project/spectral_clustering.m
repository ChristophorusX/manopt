function [ clustering, eigenvalues, eigenvectors ] = spectral_clustering( G )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform spectral clustering function on the graph, trying to recover the
% planted partition by grouping the eigenvectors of the adjacency matrix
% with k-means. (The helping function is in seperate function file.)
% @author: Ruitu Xu
%
% @parameter: A matrix G to perform spectral clustering on.
%
% @return: The resulting clustering of all the vertices in -1 and +1.
% @return: All the eigenvalues of the matrix.
% @return: All the eigenvectors of the matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[groups,eigenvalues,eigenvectors] = SpectralClustering(G,2); % Seperate into 2 groups
clustering = 2*(groups(:,2)-1)-1;
% display(clustering)
disp('Spectral clustering finished!')

end
