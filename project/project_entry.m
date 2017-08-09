%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and Z2 synchronization problem. It also
% evaluates the performance of Burer-Monteiro approach under different
% models.
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the SBM model to use
n = 5000;
a = 28;
b = 12;
% Generate adjacency matrix
[A, z_sbm] = generate_sbm_adjacency(n, a, b);
% display(A)
% Put matrix A into monotone adversary
[ V, good_vector, tagged_vector ] = monotone_adversary( A, a, b );
% display(sum(tagged_vector))
% display(sum(sum(A-V)))
% display(good_vector)
% Perform Burer-Monteiro on A and V seperately
warning('off', 'manopt:getHessian:approx')
warning('off', 'manopt:elliptopefactory:exp')
[ Q_A, Q_Acost, info_A, options_A ] = burer_monteiro( A );
[ Q_V, Q_Vcost, info_V, options_V ] = burer_monteiro( V );
[ true_cost_value_A, correlation_A ] = evaluate_performance( z_sbm, A, Q_A );
[ true_cost_value_V, correlation_V ] = evaluate_performance( z_sbm, V, Q_V );
disp(['The output optimal cost by Burer-Monteiro on A is ' num2str(Q_Acost) '.'])
disp(['The planted cost by Burer-Monteiro on A is ' num2str(true_cost_value_A) '.'])
disp(['The correlation between output X=Q_AQ_A^T and the planted vector z is ' ...
    num2str(correlation_A) '.'])
disp(['The output optimal cost by Burer-Monteiro is ' num2str(Q_Vcost) '.'])
disp(['The planted cost by Burer-Monteiro on V is ' num2str(true_cost_value_V) '.'])
disp(['The correlation between output X=Q_VQ_V^T and the planted vector z is ' ...
    num2str(correlation_V) '.'])
figure;
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
