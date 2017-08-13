%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and trying to break the sdp method
% and Burer-Monteiro approach on it.
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1000;
a = 15;
b = 2;
rng('shuffle');
[ A, z_sbm ] = generate_sbm_adjacency_logrithmic( n, a, b );
all_ones = ones(n, 1);
B = 2 * A - (all_ones * all_ones' + eye(n));
[ X, Xcost, D, eigenvalues ] = sdp_solver_sym(B);
disp(['Check dual optimal: ' num2str(trace(D))]);
D_planted_pre = diag(z_sbm) * B * diag(z_sbm);
D_planted = diag(sum(D_planted_pre, 2));
e = eig(D_planted - B);
eigens = sort(e);
% display(eigenvalues)
disp('The smallest 10 eigenvalues are: ');
for iter = 1:10
    disp(num2str(eigens(iter)));
end
