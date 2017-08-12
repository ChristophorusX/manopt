%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and trying to break the spectral method
% on it (in sparse regime).
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 2000;
% [ ~, z_sbm ] = generate_sbm_adjacency(n,a,b);
% Fix a percentage of partial recovery.
percentage = 0.5;
% Set a to a fixed number and the set b according to the treshold.
a = 25;
% Set an array recording the times of failing trails.
num_of_trails = 15;
num_of_repititions = 100;
density_of_jump = 1;
fail_records_A = zeros(num_of_trails,1);
fail_records_V_MPW = zeros(num_of_trails,1);
fail_records_V_hub = zeros(num_of_trails,1);
fail_records_V_random = zeros(num_of_trails,1);
for iteration = 1:num_of_trails
    disp(['Working on iteration ' num2str(iteration) '.'])
    jump = iteration/density_of_jump;
    % obtain b from a and threshold.
    b = ((2*a+(2+jump))-sqrt((2*a+(2+jump))^2-4*a*(a-(2+jump))))/2;
    disp(['b = ' num2str(b)])
    % A_store = zeros(n,n,num_of_repititions);
    % V_MPW_store = zeros(n,n,num_of_repititions);
    % for iter = 1:num_of_repititions
    %     A_store(:,:,iter) = generate_sbm_adjacency(n,a,b);
    %     V_MPW_store(:,:,iter) = monotone_adversary_MPW(A_store(:,:,iter),a,b);
    % end
    % Generate adjacency matrix.
    for repitition = 1:num_of_repititions
        rng('shuffle');
        disp(['Working on repitition ' num2str(repitition) '.'])
        [ A, z_sbm ] = generate_sbm_adjacency(n, a, b);
        % Put matrix A into monotone adversary.
        [ V_MPW, good_vector, tagged_vector ] = monotone_adversary_MPW( A, a, b );
        V_hub = monotone_adversary_hub( A, delta );
        V_random = monotone_adversary_random( A, delta );
        % display(sum(tagged_vector))
        % display(sum(sum(A-V_MPW)))
        % display(sum(good_vector))
        % [ clustering_A, eigenvalues_A, eigenvectors_A ] = ...
        % spectral_clustering( A_store(:,:,repitition) );
        % [ clustering_V_MPW, eigenvalues_V_MPW, eigenvectors_V_MPW ] = ...
        % spectral_clustering( V_MPW_store(:,:,repitition) );
        [ clustering_A, eigenvalues_A, eigenvectors_A ] = spectral_clustering( A );
        [ clustering_V_MPW, eigenvalues_V_MPW, eigenvectors_V_MPW ] = spectral_clustering( V_MPW );
        [ clustering_V_hub, eigenvalues_V_hub, eigenvectors_V_hub ] = spectral_clustering( V_hub );
        [ clustering_V_random, eigenvalues_V_random, eigenvectors_V_random ] = spectral_clustering( V_random );
        error_rate_A = compute_error_rate(clustering_A,z_sbm);
        error_rate_V_MPW = compute_error_rate(clustering_V_MPW,z_sbm);
        error_rate_V_hub = compute_error_rate(clustering_V_hub,z_sbm);
        error_rate_V_random = compute_error_rate(clustering_V_random,z_sbm);
        disp(['Error rate for recovery from A is ' num2str(error_rate_A)])
        disp(['Error rate for recovery from V_MPW is ' num2str(error_rate_V_MPW)])
        disp(['Error rate for recovery from V_hub is ' num2str(error_rate_V_hub)])
        disp(['Error rate for recovery from V_random is ' num2str(error_rate_V_random)])
        if error_rate_A < (1 - percentage)
            disp(['GOOD RESULT on random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_A(iteration) = fail_records_A(iteration) + 1;
        end
        if error_rate_V_MPW < (1 - percentage)
            disp(['GOOD RESULT on semi-random model V_MPW with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on semi-random model V_MPW with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_V_MPW(iteration) = fail_records_V_MPW(iteration) + 1;
        end
        if error_rate_V_hub < (1 - percentage)
            disp(['GOOD RESULT on semi-random model V_hub with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on semi-random model V_hub with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_V_hub(iteration) = fail_records_V_hub(iteration) + 1;
        end
        if error_rate_V_random < (1 - percentage)
            disp(['GOOD RESULT on semi-random model V_random with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on semi-random model V_random with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_V_random(iteration) = fail_records_V_random(iteration) + 1;
        end
        disp(['Finishing on repitition ' num2str(repitition) '.'])
    end
    disp(['Finishing on iteration ' num2str(iteration) '.'])
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Summary');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
for iteration = 1:num_of_trails
    jump = iteration/density_of_jump;
    disp(['The trail failed ' num2str(fail_records_A(iteration)) ' times ' ...
    'on random model with gap ' num2str(jump) ' from threshold.']);
    disp(['The trail failed ' num2str(fail_records_V_MPW(iteration)) ' times ' ...
    'on semi-random model V_MPW with gap ' num2str(jump) ' from threshold.']);
    disp(['The trail failed ' num2str(fail_records_V_hub(iteration)) ' times ' ...
    'on semi-random model V_hub with gap ' num2str(jump) ' from threshold.']);
    disp(['The trail failed ' num2str(fail_records_V_random(iteration)) ' times ' ...
    'on semi-random model V_random with gap ' num2str(jump) ' from threshold.']);
end
