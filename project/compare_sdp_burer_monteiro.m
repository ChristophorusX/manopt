%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and Z2 synchronization problem.
% We are trying to break the sdp method and Burer-Monteiro approach on
% these problems.
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'manopt:getHessian:approx')
warning('off', 'manopt:elliptopefactory:exp')
n = 2000;
a = 25;
num_of_trails = 10;
num_of_repititions = 10;
density_of_jump = 1;
lambda_base = 10;
percent_of_elements_being_one = 0.5;
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
disp('%%%%% Starting comparing on synchronization model with Gaussian noise. %%%%%');
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
for iter = 1:num_of_trails
    jump = iter/density_of_jump;
    lambda = lambda_base + jump;
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    disp(['>>>>>>>>This is trail ' num2str(iter) ' with lambda ' num2str(lambda)]);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    for subiter = 1:num_of_repititions
        [ Y, z_syn ] = generate_synchronization_gaussian( n, percent_of_elements_being_one, lambda );
        % Perform Burer-Monteiro on Y.
        [ Q_Y, Q_Ycost, info_Y, options_Y ] = burer_monteiro( Y );
        % Evaluate the performance of Burer-Monteiro on model.
        [ true_cost_value_Y, correlation_Y ] = evaluate_performance( z_syn, Y, Q_Y );
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Start printing evaluation output on BM-Synchronization.')
        disp('=========================================================================')
        disp(['The output optimal cost by Burer-Monteiro on Y is ' num2str(-Q_Ycost) '.'])
        disp(['The planted cost by Burer-Monteiro on Y is ' num2str(true_cost_value_Y) '.'])
        disp(['The correlation between output X=Q_YQ_Y^T and the planted vector z is ' ...
            num2str(correlation_Y) '.'])
        clustering_Y = k_means_rows(Q_Y);
        error_rate_Y = compute_error_rate(clustering_Y,z_syn);
        disp(['Error rate for recovery from Y is ' num2str(error_rate_Y)])

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Exact recovery??????')
        if correlation_Y > eps
            disp('Yes, it achieves exact recovery.');
        else
            disp('No, it does not achieve exact recovery.');
        end
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        % [ ~, Xcost, ~, eigenvalues, D_dual ] = sdp_solver_sym( Y, z_syn );
        % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        % disp('Start printing evaluation output on SDP-Synchronization.')
        % disp('=========================================================================')
        % disp(['The output primal optimum is: ' num2str(Xcost)]);
        % disp(['Check dual optimality: ' num2str(trace(D_dual))]);
        % % display(eigenvalues)
        % disp('The smallest 10 eigenvalues are: ');
        % for subsubiter = 1:10
        %     disp(num2str(eigenvalues(subsubiter)));
        % end
        % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        % disp('Exact recovery??????')
        % if eigenvalues(2) > eps
        %     disp('Yes, it achieves exact recovery.');
        % else
        %     disp('No, it does not achieve exact recovery.');
        % end
        % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
end


% b_base = 10;
% disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
% disp('%%%%%%%%%%%%%%%%%% Starting comparing on SBM model. %%%%%%%%%%%%%%%%%%');
% disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
% for iter = 1:num_of_trails
%     jump = iter/density_of_jump;
%     b = b_base + jump;
%     disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
%     disp(['>>>>>>>>This is trail ' num2str(iter) ' with a ' num2str(a) ' and b ' num2str(b)]);
%     disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
%     for subiter = 1:num_of_repititions
%         rng('shuffle');
%         [ A, z_sbm ] = generate_sbm_adjacency_logrithmic( n, a, b );
%         [ Q_A, Q_Acost, info_A, options_A ] = burer_monteiro( A );
%         [ true_cost_value_A, correlation_A ] = evaluate_performance( z_sbm, A, Q_A );
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         disp('Start printing evaluation output on BM-SBM.')
%         disp('=========================================================================')
%         disp(['The output optimal cost by Burer-Monteiro on A is ' num2str(-Q_Acost) '.'])
%         disp(['The planted cost by Burer-Monteiro on A is ' num2str(true_cost_value_A) '.'])
%         disp(['The correlation between output X=Q_AQ_A^T and the planted vector z is ' ...
%             num2str(correlation_A) '.'])
%         clustering_Y = k_means_rows(Q_Y);
%         error_rate_Y = compute_error_rate(clustering_Y,z_sbm);
%         disp(['Error rate for recovery from Y is ' num2str(error_rate_Y)])
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         disp('Exact recovery??????')
%         if correlation_A > eps
%             disp('Yes, it achieves exact recovery.');
%         else
%             disp('No, it does not achieve exact recovery.');
%         end
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%
%         all_ones = ones(n, 1);
%         B = 2 * A - (all_ones * all_ones' + eye(n));
%         [ ~, Xcost, ~, eigenvalues, D_dual ] = sdp_solver_sym(B, z_sbm);
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         disp('Start printing evaluation output on SDP-SBM.')
%         disp('=========================================================================')
%         disp(['The output primal optimum is: ' num2str(Xcost)]);
%         disp(['Check dual optimality: ' num2str(trace(D_dual))]);
%         % display(eigenvalues)
%         disp('The smallest 10 eigenvalues are: ');
%         for subsubiter = 1:10
%             disp(num2str(eigenvalues(subsubiter)));
%         end
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         disp('Exact recovery??????')
%         if eigenvalues(2) > eps
%             disp('Yes, it achieves exact recovery.');
%         else
%             disp('No, it does not achieve exact recovery.');
%         end
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%     end
% end
