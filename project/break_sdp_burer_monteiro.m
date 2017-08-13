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
A = generate_sbm_adjacency_logrithmic( n, a, b );
all_ones = ones(n, 1);
B = 2 * A - (all_ones * all_ones' - eye(n));
[ X, Xcost, D, eigenvalues ] = sdp_solver_sym(B);
disp(['Check dual optimal: ' num2str(trace(D))]);
disp('The smallest 10 eigenvalues are: ');
for iter = 1:10
    disp(num2str(eigenvalues(iter)));
end
