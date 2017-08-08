%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and Z2 synchronization problem. 
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random adjacency matrix (symmetric) A of dimension n (even) 
% according to inner- and inter-cluster edge probability p and q.
% NOTE: this setting is for stochastic block model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_of_vertices = 5000;
a = 16;
b = 4;
n = number_of_vertices;
p = a*log(n)/n;
q = b*log(n)/n;
% Properties of the planted partition.
z_sbm = [ones(1,n/2) -ones(1,n/2)]'; % planted partition for SBM.
R1_pre = rand(n/2,n/2) > (1-p);
% display(R1_pre)
R3_pre = rand(n/2,n/2) > (1-p);
% display(R3_pre)
for i = 1:(n/2)
    for j = 1:i
        if j < i
            R1_pre(i,j) = R1_pre(j,i);
            R3_pre(i,j) = R3_pre(j,i);
        else
            R1_pre(i,i) = 1;
            R3_pre(i,i) = 1;
        end
    end
end
% display(R1_pre)
% display(R3_pre)
R2 = rand(n/2,n/2) > (1-q);
A = [R1_pre R2; R2' R3_pre];
% display(A)
row_sum = sum(A,2);
for i = 1:n
    disp(['The summation of row #' num2str(i) ' is ' num2str(row_sum(i)) '.'])
end
disp(['The average edge of every row should be around ' num2str((p*n + q*n)/2) '.'])
disp(['The empirical average number of edges is ' num2str(sum(row_sum)/n)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement a monotone adversary described in MPW16. It removes with
% probability delta the edges of node with only two neighbors of different
% spin. Delta is given by epsilon, the rate of b/(a+b).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign a matrix A as the adjacency matrix for the adversary.
with_diagonal = 0;
diagonal_elements = diag(A);
if diagonal_elements ~= 0
    V = A-diag(A); % If input A has diag(A)=1, remove the diagonal.
    with_diagonal = 1;
end
epsilon = b/(a+b);
if epsilon <= 1/3
    delta = 1;
else
    delta = (1-2*epsilon)^2/epsilon^2;
end
degree_vector = sum(V,2);
% Step 1: mark all the 'good' nodes.
good_vector = zeros(n,1);
for i = 1:n
    if degree_vector(i) >= 3
        neighbor_degree_vector = degree_vector.*V(i,:);
        degree_check_vector = (neighbor_degree_vector ~= 2).* ...
            (neighbor_degree_vector ~= 0); % Perform a filter on a filter.
        if sum(degree_check_vector) >= 3 % If at least exists 3 such neighbors.
            good_vector(i) = 1;
        end
    end
end
% Step 2: mark all the 'tagged' nodes.
tagged_vector = zeros(n,1);
for i = 1:n
    if degree_vector(i) == 2
        align_check_vector = good_vector.*V(i,:);
        if sum(align_check_vector) == 2 % See if both neighbors are good.
            tagged_vector(i) = 1;
        end
    end
end
% Step 3: remove edges of tagged nodes with probability delta.
for i = 1:n
    if tagged_vector(i) == 1
        if rand(1,1) < delta
            V(i,:) = zeros(1,n);
        end
    end
end
if with_diagonal == 1
    V = V + diag(diagonal_elements);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random matrix Y of dimension n for Z2 synchronization problem,
% where Y=zz^T+\sigma W, with W a Wigner matrix generated from Gaussian
% distribution.
% NOTE: this setting is for synchronization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
percent_of_elements_being_one = 0.5;
z_pre = rand(n,1) > (1-percent_of_elements_being_one);
z_syn = (2*z_pre)-1;
% display(z_syn)
R = normrnd(0,1,n,n);
W = R - diag(R);
% display(W)
lambda = 20; % lambda is the value to move around.
sigma = sqrt(n)/lambda;
Y = z_syn*z_syn'+ sigma*W;
% display(Y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform a manifold optimization on the cost function Tr(Q^TBQ) given the
% manifold (Symmetric positive semidefinite, fixed-rank with unit diagonal)
% with rank(X=QQ^T)<=2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = Y; % Assign a specific problem to the solver
k = 2;
manifold = elliptopefactory(n, k);
problem.M = manifold;
% Note that Manopt only minimize the cost function.
problem.cost = @(Q) -trace(Q'*B*Q);
% The Euclidian gradient is given by calculation.
problem.egrad = @(Q) -((B+B')*Q+[diag(B*diag(Q(:,1))) diag(B*diag(Q(:,2)))]);
% Numerically check gradient and Hessian consistency.
figure;
checkgradient(problem);
figure;
checkhessian(problem);
% Pass the problem to the solver, i.e. trust regions.
[Q, Qcost, info, options] = trustregions(problem);
disp(['The output optimal cost by the algorithm is ' num2str(Qcost) '.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing some evaluations on the above methods with different problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = z_syn; % Assign a specific problem to the evaluator.
% disp(z)
true_cost_value = -z'*B*z; % Cost given ground truth z.
% disp(true_cost_value)
disp(['The planted cost value of this problem is ' num2str(true_cost_value) '.'])
% Compute the correlation between the critical point Q and ground truth z.
correlation = sqrt(z'*(Q*Q')*z)/n;
disp(['The correlation between output X=QQ^T and the planted vector z is ' ...
    num2str(correlation) '.'])
% Figure on convergence of the algorithm.
figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration #');
ylabel('Gradient norm');
title('Convergence of the trust-regions algorithm on the manifold');

