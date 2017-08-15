warning('off', 'manopt:getHessian:approx')
warning('off', 'manopt:elliptopefactory:exp')
n = 1000;
a = 50;
num_of_trails = 10;
num_of_repititions = 1;
density_of_jump = 1;
lambda_base = 10;
delta = 1/10;
percent_of_elements_being_one = 0.5;
b_base = 2;
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
disp('%%%%%%%%%%%%%%%%%% Starting comparing on SBM model. %%%%%%%%%%%%%%%%%%');
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
for iter = 1:num_of_trails
    jump = iter/density_of_jump;
    b = b_base + jump;
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    disp(['>>>>>>>>This is trail ' num2str(iter) ' with a ' num2str(a) ' and b ' num2str(b)]);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    for subiter = 1:num_of_repititions
        rng('shuffle');
        [ A, z_sbm ] = generate_sbm_adjacency_logrithmic( n, a, b );
        [ V_MPW, good_vector, tagged_vector ] = monotone_adversary_MPW( A, a, b );
        V_hub = monotone_adversary_hub( A, delta );
        V_random = monotone_adversary_random( A, delta );

        [ Q_A, Q_Acost, info_A, options_A ] = burer_monteiro( A );
        [ Q_V_MPW, Q_V_MPWcost, info_V_MPW, options_V_MPW ] = burer_monteiro( V_MPW );
        [ Q_V_hub, Q_V_hubcost, info_V_hub, options_V_hub ] = burer_monteiro( V_hub );
        [ Q_V_random, Q_V_randomcost, info_V_random, options_V_random ] = burer_monteiro( V_random );

        [ true_cost_value_A, correlation_A ] = evaluate_performance( z_sbm, A, Q_A );
        [ true_cost_value_V_MPW, correlation_V_MPW ] = evaluate_performance( z_sbm, V_MPW, Q_V_MPW );
        [ true_cost_value_V_hub, correlation_V_hub ] = evaluate_performance( z_sbm, V_hub, Q_V_hub );
        [ true_cost_value_V_random, correlation_V_random ] = evaluate_performance( z_sbm, V_random, Q_V_random );
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Start printing evaluation output on BM-SBM.')
        disp('=========================================================================')
        disp(['The output optimal cost by Burer-Monteiro on A is ' num2str(-Q_Acost) '.'])
        disp(['The planted cost by Burer-Monteiro on A is ' num2str(true_cost_value_A) '.'])
        disp(['The correlation between output X=Q_AQ_A^T and the planted vector z is ' ...
            num2str(correlation_A) '.'])
        disp('=========================================================================')
        disp(['The output optimal cost by Burer-Monteiro on V_MPW is ' num2str(-Q_V_MPWcost) '.'])
        disp(['The planted cost by Burer-Monteiro on V_MPW is ' num2str(true_cost_value_V_MPW) '.'])
        disp(['The correlation between output X=Q_V_MPWQ_V_MPW^T and the planted vector z is ' ...
            num2str(correlation_V_MPW) '.'])
        disp('=========================================================================')
        disp(['The output optimal cost by Burer-Monteiro on V_hub is ' num2str(-Q_V_hubcost) '.'])
        disp(['The planted cost by Burer-Monteiro on V_hub is ' num2str(true_cost_value_V_hub) '.'])
        disp(['The correlation between output X=Q_V_hubQ_V_hub^T and the planted vector z is ' ...
            num2str(correlation_V_hub) '.'])
        disp('=========================================================================')
        disp(['The output optimal cost by Burer-Monteiro on V_random is ' num2str(-Q_V_randomcost) '.'])
        disp(['The planted cost by Burer-Monteiro on V_random is ' num2str(true_cost_value_V_random) '.'])
        disp(['The correlation between output X=Q_V_randomQ_V_random^T and the planted vector z is ' ...
            num2str(correlation_V_random) '.'])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Exact recovery??????')
        if correlation_A > eps
            disp('Yes, it achieves exact recovery.');
        else
            disp('No, it does not achieve exact recovery.');
        end
        all_ones = ones(n, 1);
        B = 2 * A - (all_ones * all_ones' + eye(n));
        D_planted_pre = diag(z_sbm) * B * diag(z_sbm);
        D_planted = diag(sum(D_planted_pre, 2));
        e = eig(D_planted - diag(z_sbm) * B * diag(z_sbm));
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
