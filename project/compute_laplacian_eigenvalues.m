function [ eigenvalues ] = compute_laplacian_eigenvalues( Y, z )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gives a list of eigenvalues of the Laplacian matrix
% D_{diag(z)Y diag(z)} - diag(z)Y diag(z).
% @author: Ruitu Xu
%
% @parameter: Y be the symmetric matrix given by the model.
% @parameter: z be the planted value of the vertices.
%
% @return: eigenvalues of the Laplacian matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_planted_pre = diag(z) * Y * diag(z);
D_planted = diag(sum(D_planted_pre, 2));
e = eig(D_planted - diag(z) * Y * diag(z));
eigenvalues_pre = sort(e);
eigenvalues = eigenvalues_pre .* (eigenvalues_pre > eps);

end
