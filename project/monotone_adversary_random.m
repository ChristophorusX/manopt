function [ V ] = monotone_adversary_random( A, delta )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement a monotone adversary creating random modify on graph.
%
% @parameter: An adjacency matrix given by a SBM model.
% @parameter: constant delta for probability of creating and deleting
% edges.
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
for iter = 1:dimension
    adding_sign = - (2 * V(iter,:) - 1)';
    first_half = [ones(1,dimension/2) zeros(1,dimension/2)]';
    second_half = [zeros(1,dimension/2) ones(1,dimension/2)]';
    if iter <= dimension/2
        row_first_half = first_half .* V(iter,:)';
        row_second_half = second_half .* V(iter,:)';
        validity_vec = ((1 - (2 * row_first_half - 1)) / 2) .* first_half + ...
        row_second_half;
    else
        row_first_half = first_half .* V(iter,:)';
        row_second_half = second_half .* V(iter,:)';
        validity_vec = ((1 - (2 * row_second_half - 1)) / 2) .* second_half + ...
        row_first_half;
    end
    rng('shuffle');
    moving_vec = rand(dimension,1) > (1 - delta);
    V(iter,:) = V(iter,:) + (moving_vec .* adding_sign .* validity_vec)';
    % for subiter = iter:dimension
    %     if V(iter,subiter) == 0 && ((iter <= dimension/2 && subiter > dimension/2) ...
    %         || (iter > dimension/2 && subiter <= dimension/2))
    %         rng('shuffle');
    %         if rand(1,1) < (1 - delta)
    %             V(iter,subiter) = 1;
    %             V(subiter,iter) = 1;
    %         end
    %     elseif V(iter,subiter) == 1 && ~((iter <= dimension/2 && subiter > dimension/2) ...
    %         || (iter > dimension/2 && subiter <= dimension/2))
    %         rng('shuffle');
    %         if rand(1,1) < delta
    %             V(iter,subiter) = 0;
    %             V(subiter,iter) = 0;
    %         end
    %     end
    % end
end
for iter = 1:dimension
    for subiter = iter:dimension
        V(subiter,iter) = V(iter,subiter);
    end
end
if with_diagonal == 1
    V = V - diag(diag(V)) + diag(diagonal_elements);
end
% display(V)
disp('Monotone adversary (random version) performed!')

end
