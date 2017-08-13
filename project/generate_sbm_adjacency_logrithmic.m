function [ A, z_sbm ] = generate_sbm_adjacency_logrithmic( number_of_vertices, a, b )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random adjacency matrix (symmetric) A of dimension n (even)
% according to inner- and inter-cluster edge probability p and q.
% NOTE: this setting is for stochastic block model in logrithmic regime.
%
% @parameter: number of vertices in the model.
% @parameter: constant a for inner cluster probability.
% @parameter: constant b for inter cluster probability.
%
% @return: The adjacency matrix A corresponding to the random graph
% generated by the SBM of two clusters.
% @return: The ground truth z_sbm marked by +1 and -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = number_of_vertices;
p = a*log(n)/n;
q = b*log(n)/n;
% Properties of the planted partition.
z_sbm = [ones(1,n/2) -ones(1,n/2)]'; % planted partition for SBM.
R1_pre = rand(n/2,n/2) > (1-p);
% display(R1_pre)
R3_pre = rand(n/2,n/2) > (1-p);
% display(R3_pre)
R1 = zeros(n/2);
R3 = zeros(n/2);
for i = 1:(n/2)
    for j = 1:(n/2)
        if j > i
            R1(i,j) = R1_pre(i,j);
            R3(i,j) = R3_pre(i,j);
        elseif j < i
            R1(i,j) = R1_pre(j,i);
            R3(i,j) = R3_pre(j,i);
        else
            R1(i,i) = 1;
            R3(i,i) = 1;
        end
    end
end
% display(R1_pre)
% display(R3_pre)
R2 = rand(n/2,n/2) > (1-q);
A = [R1 R2; R2' R3]; % <= NOTE: A is the adjacency matrix of SBM
% display(A)
% row_sum = sum(A,2);
% for i = 1:n
%     disp(['The summation of row #' num2str(i) ' is ' num2str(row_sum(i)) '.'])
% end
% disp(['The average edge of every row should be around ' num2str((p*n + q*n)/2) '.'])
% disp(['The empirical average number of edges is ' num2str(sum(row_sum)/n)])
disp('Adjacency matrix for SBM generated!')

end
