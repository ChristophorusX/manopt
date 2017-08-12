function [ V, good_vector, tagged_vector ] = monotone_adversary( A, a, b )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement a monotone adversary described in MPW16. It removes with
% probability delta the edges of node with only two neighbors of different
% spin. Delta is given by epsilon, the rate of b/(a+b).
%
% @parameter: An adjacency matrix given by a SBM model.
% @parameter: constant a for inner cluster probability.
% @parameter: constant b for inter cluster probability.
%
% @return: An adjacency matrix V after the operation of the adversary.
% @return: An 0-1 vector giving all the vertices marked 'good'.
% @return: An 0-1 vector giving all the vertices marked 'tagged'.
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
epsilon = b/(a+b);
if epsilon <= 1/3
    delta = 1;
else
    delta = (1-2*epsilon)^2/epsilon^2;
end
% display(delta)
degree_vector = sum(V,2);
% display(degree_vector)
% Step 1: mark all the 'good' nodes.
[ dimension, ~ ] = size(A);
good_vector = zeros(dimension,1);
for iter = 1:dimension
    if degree_vector(iter) >= 3
        neighbor_degree_vector = degree_vector.*V(iter,:)';
        degree_check_vector = (neighbor_degree_vector ~= 2).* ...
            (neighbor_degree_vector ~= 0); % Perform a filter on a filter.
        if sum(degree_check_vector) >= 3 % If at least exists 3 such neighbors.
            good_vector(iter) = 1;
        end
    end
end
% Step 2: mark all the 'tagged' nodes.
tagged_vector = zeros(dimension,1);
for iter = 1:dimension
    if degree_vector(iter) == 2
        align_check_vector = good_vector.*V(iter,:)';
        if sum(align_check_vector) == 2 % See if both neighbors are good.
            tagged_vector(iter) = 1;
        end
    end
end
% Step 3: remove edges of tagged nodes with probability delta.
for iter = 1:dimension
    if tagged_vector(iter) == 1
        rng('shuffle');
        if rand(1,1) < delta
            V(iter,:) = V(iter,:) - V(iter,:);
            V(:,iter) = V(:,iter) - V(iter,:)';
        end
    end
end
if with_diagonal == 1
    V = V + diag(diagonal_elements); % <= NOTE: V is adjacency matrix after adversary
end
% display(V)
disp('Monotone adversary performed!')

end
