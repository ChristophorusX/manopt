% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters. @author: Ruitu Xu

% Generate random adjacency matrix (symmetric) A of dimension n (even) 
% according to inner- and inter-cluster edge probability p and q

number_of_vertices = 5000;
a = 16;
b = 4;
n = number_of_vertices;
p = a*log(n)/n;
q = b*log(n)/n;
% Properties of the planted partition
z_sbm = [ones(1,n/2) -ones(1,n/2)]'; % planted partition
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

% Generate random matrix Y of dimension n for Z2 synchronization problem,
% where Y=zz^T+\sigma W, with W a Wigner matrix generated from Gaussian
% distribution

percent_of_elements_being_one = 0.5;
z_pre = rand(n,1) > (1-percent_of_elements_being_one);
z_syn = (2*z_pre)-1;
% display(z_syn)

R = normrnd(0,1,n,n);
W = R - diag(R);
% display(W)

lambda = 20; % value to move around
sigma = sqrt(n)/lambda;
Y = z*z'+ sigma*W;
% display(Y)

% Perform a manifold optimization on the cost function Tr(Q^TBQ) given the
% manifold (Symmetric positive semidefinite, fixed-rank with unit diagonal)
% with rank(X=QQ^T)<=2

B = Y; % Assign a specific problem to the solver
k = 2;
manifold = elliptopefactory(n, k);
problem.M = manifold;

problem.cost = @(Q) -trace(Q'*B*Q); % Note that Manopt only minimize the cost function
problem.egrad = @(Q) -((B+B')*Q+[diag(B*diag(Q(:,1))) diag(B*diag(Q(:,2)))]); % by calculation

% Numerically check gradient and Hessian consistency
figure;
checkgradient(problem);
figure;
checkhessian(problem);

[Q, Qcost, info, options] = trustregions(problem);
disp(['The output optimal cost by the algorithm is ' num2str(Qcost) '.'])

z = z_syn; % Assign a specific problem to the evaluator
% disp(z)
true_cost_value = -z'*B*z;
% disp(true_cost_value)
disp(['The planted cost value of this problem is ' num2str(true_cost_value) '.'])

correlation = sqrt(z'*(Q*Q')*z)/n;
disp(['The correlation between output X=QQ^T and the planted vector z is ' num2str(correlation) '.'])

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration #');
ylabel('Gradient norm');
title('Convergence of the trust-regions algorithm on the manifold');

