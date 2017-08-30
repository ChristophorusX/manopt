import matlab.engine
import numpy as np
import scipy as sp
ml = matlab.engine.start_matlab()

ml.warning('off', 'manopt:getHessian:approx', nargout=0)
ml.warning('off', 'manopt:elliptopefactory:exp', nargout=0)
n = 1000
n_ml = ml.double(n)
num_of_trails = 10
num_of_repititions = 1
density_of_jump = 10
lambda_base = np.sqrt(2 * np.log(n))
delta = 1 / 10
percent_of_elements_being_one = 0.5
percent_of_elements_being_one_ml = ml.double(percent_of_elements_being_one)
print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
print('%%%%%%%%%%%%%%%% Starting comparing on sync model. %%%%%%%%%%%%%%')
print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
for iter in range(0, num_of_trails):
    jump = iter / density_of_jump
    lambda_SNR = lambda_base * (3 + jump)
    lambda_SNR_ml = ml.double(float(lambda_SNR))
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    print('>>>>>>>>This is trail' + str(iter) +
          ' with lambda ' + str(lambda_SNR))
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    for subiter in range(0, num_of_repititions):
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print('>>>>>>>>Rep #' + str(subiter) +
              ' with lambda ' + str(lambda_SNR))
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        ml.rng('shuffle', nargout=0)
        Y_normalized, z_syn = ml.generate_synchronization_gaussian_normalized(n_ml, percent_of_elements_being_one_ml, lambda_SNR_ml, nargout=2)
        Q_Y_normalized, Q_Y_normalizedcost, info_Y_normalized, options_Y_normalized = ml.burer_monteiro(Y_normalized, nargout=4)
        true_cost_value_Y_normalized, correlation_Y_normalized = ml.evaluate_performance(z_syn, Y_normalized, Q_Y_normalized, nargout=3)
        clustering_Y_normalized = ml.k_means_rows(Q_Y_normalized)
        error_rate_Y_normalized = ml.compute_error_rate(clustering_Y_normalized, z_syn)
        print('Error rate for recovery from Y_normalized is ' +
              str(error_rate_Y_normalized))

        # plot_x_Y_normalized = Q_Y_normalized(:,1)
        # plot_y_Y_normalized = Q_Y_normalized(:,2)
        # plot(plot_x_Y_normalized,plot_y_Y_normalized, 'o', 'color', 'blue')
        # hold on
        # xlabel('x')
        # ylabel('y')
        # title('Row vectors of second order critical point \
        #  Q_Y_normalized on unit circle.')
        # hold off

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('Start printing evaluation output on BM-sync.')
        print('=======================================================')
        print('The output optimal cost by Burer-Monteiro on Y_normalized is ' +
              str(-Q_Y_normalizedcost) + '.')
        print('The planted cost by Burer-Monteiro on Y_normalized is ' +
              str(true_cost_value_Y_normalized) + '.')
        print('The correlation between output X=Q_Y_normalizedQ_Y_normalized^T \
              and the planted vector z is ' +
              str(correlation_Y_normalized) + '.')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('Exact recovery??????')
        if correlation_Y_normalized > eps:
            print('Yes, it achieves exact recovery.')
        else:
            print('No, it does not achieve exact recovery.')
        # all_ones = ones(n, 1)
        # B = 2 * Y_normalized - (all_ones * all_ones' + eye(n))
        # D_planted_pre = diag(z_syn) * B * diag(z_syn)
        # D_planted = diag(sum(D_planted_pre, 2))
        # e = eig(D_planted - diag(z_syn) * B * diag(z_syn))
        # eigenvalues_pre = sort(e)
        # eigenvalues = eigenvalues_pre .* (eigenvalues_pre > eps)
        # print('The smallest 10 eigenvalues are: ')
        # for subsubiter in range(1, 10):
        #     print(eigenvalues(subsubiter))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
