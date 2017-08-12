function [ V ] = monotone_adversary_hub( A, delta )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement a monotone adversary creating a hub on graph. Add to the
% first point a linear portion of the total possible edges within the
% group and remove a linar portion of those crossing the groups.
%
% @parameter: An adjacency matrix given by a SBM model.
% @parameter: constant delta for probability of creating hub edges.
%
% @return: An adjacency matrix V after the operation of the adversary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign a matrix A as the adjacency matrix for the adversary.
with_diagonal = 0;
diagonal_elements = diag(A);
if diagonal_elements ~= 0
    V = A-diag(diag(A)); % If input A has diag(A)=1, remove the diagonal.
    with_diagonal = 1;
else
    V = A;
end
% display(V)
% display(delta)
[ dimension, ~ ] = size(A);
for iter = 2:dimension
    if V(iter,1) == 0 && iter <= dimension/2
        rng('shuffle');
        if rand(1,1) < (1 - delta)
            V(iter,1) = 1;
        end
    elseif V(iter,1) == 1 && iter > dimension/2
        rng('shuffle');
        if rand(1,1) < delta
            V(iter,1) = 0;
        end
    end
end
V(1,:) = V(:,1);
if with_diagonal == 1
    V = V + diag(diagonal_elements); % <= NOTE: V is adjacency matrix after adversary
end
% display(V)
disp('Monotone adversary (hub version) performed!')

end
