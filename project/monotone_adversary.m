% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters. @author: Ruitu Xu

% Generate random adjacency matrix (symmetric) of dimension n (even) 
% according to inner- and inter-cluster edge probability p and q

number_of_vertices = 30;
a = 16;
b = 4;
n = number_of_vertices;
p = a*log(n)/n;
q = b*log(n)/n;
R1_pre = rand(n/2,n/2) > (1-p);
%display(R1_pre)
R3_pre = rand(n/2,n/2) > (1-p);
%display(R3_pre)
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
%display(R1_pre)
%display(R3_pre)

R2 = rand(n/2,n/2) > (1-q);
A = [R1_pre R2; R2' R3_pre];
%display(A)
row_sum = sum(A,2);
for i = 1:n
    disp(['The summation of row #' num2str(i) ' is ' num2str(row_sum(i)) '.'])
end
disp(['The average edge of every row should be around ' num2str((p*n + q*n)/2) '.'])
disp(['The empirical average number of edges is ' num2str(sum(row_sum)/n)])

% Properties of the planted partition
z = [ones(n/2) -ones(n/2)]'; % planted partition
true_cost_value = z'*A*z;
disp(['The planted cost value of this problem is ' true_cost_value '.'])


% Perform a manifold optimization on the cost function Tr(AX) given the
% manifold (Symmetric positive semidefinite, fixed-rank with unit diagonal)
% with rank(X)<=2

k = 2;
manifold = elliptopefactory(n, k);
problem.M = manifold;

problem.cost = @(X) -sum(diag((A*X)));

[X, Xcost, info] = trustregions(problem);
disp(['The output optimal cost by the algorithm is ' Xcost '.'])

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration #');
ylabel('Gradient norm');
title('Convergence of the trust-regions algorithm on the manifold');

