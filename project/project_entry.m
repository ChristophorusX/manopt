%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and Z2 synchronization problem. It also
% evaluates the performance of Burer-Monteiro approach under different
% models.
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the SBM model to use.
n = 1000;
a = 28;
b = 12;
% Generate adjacency matrix.
[ A, z_sbm ] = generate_sbm_adjacency_sparse(n, a, b);
% display(A)
% Put matrix A into monotone adversary.
[ V, good_vector, tagged_vector ] = monotone_adversary_MPW( A, a, b );
% display(sum(tagged_vector))
% display(sum(sum(A-V)))
% display(good_vector)
% Perform Burer-Monteiro on A and V seperately.
warning('off', 'manopt:getHessian:approx')
warning('off', 'manopt:elliptopefactory:exp')
[ Q_A, Q_Acost, info_A, options_A ] = burer_monteiro( A );
[ Q_V, Q_Vcost, info_V, options_V ] = burer_monteiro( V );
% Evaluate the performance of Burer-Monteiro on model.
[ true_cost_value_A, correlation_A ] = evaluate_performance( z_sbm, A, Q_A );
[ true_cost_value_V, correlation_V ] = evaluate_performance( z_sbm, V, Q_V );
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Start printing evaluation output on BM-SBM.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['The output optimal cost by Burer-Monteiro on A is ' num2str(Q_Acost) '.'])
disp(['The planted cost by Burer-Monteiro on A is ' num2str(true_cost_value_A) '.'])
disp(['The correlation between output X=Q_AQ_A^T and the planted vector z is ' ...
    num2str(correlation_A) '.'])
disp('compared to')
disp(['The output optimal cost by Burer-Monteiro is ' num2str(Q_Vcost) '.'])
disp(['The planted cost by Burer-Monteiro on V is ' num2str(true_cost_value_V) '.'])
disp(['The correlation between output X=Q_VQ_V^T and the planted vector z is ' ...
    num2str(correlation_V) '.'])
critical_point_plot_SBM = figure('Visible','off');
plot_x_A = Q_A(:,1);
% display(plot_x_A)
plot_y_A = Q_A(:,2);
% display(plot_y_A)
plot_x_V = Q_V(:,1);
% display(plot_x_V)
plot_y_V = Q_V(:,2);
% display(plot_y_V)
plot(plot_x_A,plot_y_A, 'o', 'color', 'blue')
hold on;
plot(plot_x_V,plot_y_V, 'o', 'color', 'red')
lgd = legend('(q_i^A)^T', '(q_i^V)^T');
lgd.FontSize = 12;
xlabel('x')
ylabel('y')
title('Row vectors of second order critical point Q_A and Q_V on unit circle.')
hold off;
saveas(critical_point_plot_SBM,'SBM critical point plot','png')

% Specify the synchronization model to use.
percent_of_elements_being_one = 0.5;
lambda = 8;
[ Y, z_syn ] = generate_synchronization_matrix( n, percent_of_elements_being_one, lambda );
% Perform Burer-Monteiro on Y.
[ Q_Y, Q_Ycost, info_Y, options_Y ] = burer_monteiro( Y );
% Evaluate the performance of Burer-Monteiro on model.
[ true_cost_value_Y, correlation_Y ] = evaluate_performance( z_syn, Y, Q_Y );
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Start printing evaluation output on BM-Synchronization.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['The output optimal cost by Burer-Monteiro on Y is ' num2str(Q_Ycost) '.'])
disp(['The planted cost by Burer-Monteiro on Y is ' num2str(true_cost_value_Y) '.'])
disp(['The correlation between output X=Q_YQ_Y^T and the planted vector z is ' ...
    num2str(correlation_Y) '.'])
critical_point_plot_syn = figure('Visible','off');
plot_x_Y = Q_Y(:,1);
plot_y_Y = Q_Y(:,2);
plot(plot_x_Y,plot_y_Y, 'o', 'color', 'green')
hold on;
lgd_syn = legend('(q_i^Y)^T');
lgd_syn.FontSize = 12;
xlabel('x')
ylabel('y')
title('Row vectors of second order critical point Q_Y on unit circle.')
hold off;
saveas(critical_point_plot_syn,'Synchronization critical point plot','png')

% Apply spectral clustering on the same SBM model.
[ clustering_A, eigenvalues_A, eigenvectors_A ] = spectral_clustering( A );
[ clustering_V, eigenvalues_V, eigenvectors_V ] = spectral_clustering( V );
% Evaluate the performance of spectral method.
error_rate_A = compute_error_rate(clustering_A,z_sbm);
error_rate_V = compute_error_rate(clustering_V,z_sbm);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Start printing evaluation output on Spectral-SBM.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['The error rate of clustering on A compared to the planted vector z is ' ...
    num2str(error_rate_A) '.'])
disp(['The error rate of clustering on V compared to the planted vector z is ' ...
    num2str(error_rate_V) '.'])

% Apply spectral clustering on the same Synchronization model.
[ clustering_Y, eigenvalues_Y, eigenvectors_Y ] = spectral_clustering_rw( Y, 2 );
% Evaluate the performance of spectral method.
error_rate_Y = compute_error_rate(clustering_Y,z_syn);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Start printing evaluation output on Spectral-Synchronization.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['The error rate of clustering on Y compared to the planted vector z is ' ...
    num2str(error_rate_Y) '.'])
