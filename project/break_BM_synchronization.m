warning('off', 'manopt:getHessian:approx')
warning('off', 'manopt:elliptopefactory:exp')
n = 2000;
num_of_trails = 10;
num_of_repititions = 1;
density_of_jump = 10;
lambda_base = sqrt(2*log(n));
delta = 1/10;
percent_of_elements_being_one = 0.5;
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
disp('%%%%%%%%%%%%%%%%%% Starting comparing on sync model. %%%%%%%%%%%%%%%%%%');
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
for iter = 1:num_of_trails
    jump = iter / density_of_jump;
    lambda = lambda_base * (3 + jump);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    disp(['>>>>>>>>This is trail ' num2str(iter) ' with lambda ' num2str(lambda) ]);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    for subiter = 1:num_of_repititions
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp(['>>>>>>>>Rep #' num2str(subiter) ' with lambda ' num2str(lambda) ]);
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        rng('shuffle');
        [ Y_normalized, z_syn ] = generate_synchronization_gaussian_normalized( n, percent_of_elements_being_one, lambda );
        [ Q_Y_normalized, Q_Y_normalizedcost, info_Y_normalized, options_Y_normalized ] = burer_monteiro( Y_normalized );
        [ true_cost_value_Y_normalized, correlation_Y_normalized ] = evaluate_performance( z_syn, Y_normalized, Q_Y_normalized );
        clustering_Y_normalized = k_means_rows(Q_Y_normalized);
        error_rate_Y_normalized = compute_error_rate(clustering_Y_normalized,z_syn);
        disp(['Error rate for recovery from Y_normalized is ' num2str(error_rate_Y_normalized)])

        plot_x_Y_normalized = Q_Y_normalized(:,1);
        plot_y_Y_normalized = Q_Y_normalized(:,2);
        plot(plot_x_Y_normalized,plot_y_Y_normalized, 'o', 'color', 'blue')
        hold on;
        xlabel('x')
        ylabel('y')
        title('Row vectors of second order critical point Q_Y_normalized on unit circle.')
        hold off;

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Start printing evaluation output on BM-sync.')
        disp('=========================================================================')
        disp(['The output optimal cost by Burer-Monteiro on Y_normalized is ' num2str(-Q_Y_normalizedcost) '.'])
        disp(['The planted cost by Burer-Monteiro on Y_normalized is ' num2str(true_cost_value_Y_normalized) '.'])
        disp(['The correlation between output X=Q_Y_normalizedQ_Y_normalized^T and the planted vector z is ' ...
            num2str(correlation_Y_normalized) '.'])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Exact recovery??????')
        if correlation_Y_normalized > eps
            disp('Yes, it achieves exact recovery.');
        else
            disp('No, it does not achieve exact recovery.');
        end
        all_ones = ones(n, 1);
        B = 2 * Y_normalized - (all_ones * all_ones' + eye(n));
        D_planted_pre = diag(z_syn) * B * diag(z_syn);
        D_planted = diag(sum(D_planted_pre, 2));
        e = eig(D_planted - diag(z_syn) * B * diag(z_syn));
        eigenvalues_pre = sort(e);
        eigenvalues = eigenvalues_pre .* (eigenvalues_pre > eps);
        disp('The smallest 10 eigenvalues are: ');
        for subsubiter = 1:10
            disp(num2str(eigenvalues(subsubiter)));
        end
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
end
