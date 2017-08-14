%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and trying to break the sdp method
% and Burer-Monteiro approach on it.
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1000;
a = 15;
b = 5;
rng('shuffle');
[ A, z_sbm ] = generate_sbm_adjacency_logrithmic( n, a, b );
[ Q_A, Q_Acost, info_A, options_A ] = burer_monteiro( A );
[ true_cost_value_A, correlation_A ] = evaluate_performance( z_sbm, A, Q_A );
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Start printing evaluation output on BM-SBM.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['The output optimal cost by Burer-Monteiro on A is ' num2str(-Q_Acost) '.'])
disp(['The planted cost by Burer-Monteiro on A is ' num2str(true_cost_value_A) '.'])
disp(['The correlation between output X=Q_AQ_A^T and the planted vector z is ' ...
    num2str(correlation_A) '.'])

all_ones = ones(n, 1);
B = 2 * A - (all_ones * all_ones' + eye(n));
[ X, Xcost, D_planted, eigenvalues, D_dual ] = sdp_solver_sym(B, z_sbm);
disp(['Check dual optimality: ' num2str(trace(D_dual))]);
% display(eigenvalues)
disp('The smallest 10 eigenvalues are: ');
for iter = 1:10
    disp(num2str(eigenvalues(iter)));
end
