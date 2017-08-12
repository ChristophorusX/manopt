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
fail_records_A = zeros(num_of_trails,1);
fail_records_V = zeros(num_of_trails,1);
for iteration = 14:num_of_trails
    disp(['Working on iteration ' num2str(iteration) '.'])
    jump = iteration/1;
    % obtain b from a and threshold.
    b = ((2*a+(2+jump))-sqrt((2*a+(2+jump))^2-4*a*(a-(2+jump))))/2;
    disp(['b = ' num2str(b)])
    % A_store = zeros(n,n,num_of_repititions);
    % V_store = zeros(n,n,num_of_repititions);
    % for iter = 1:num_of_repititions
    %     A_store(:,:,iter) = generate_sbm_adjacency(n,a,b);
    %     V_store(:,:,iter) = monotone_adversary(A_store(:,:,iter),a,b);
    % end
    % Generate adjacency matrix.
    for repitition = 1:num_of_repititions
        rng('shuffle');
        disp(['Working on repitition ' num2str(repitition) '.'])
        [ A, z_sbm ] = generate_sbm_adjacency(n, a, b);
        % Put matrix A into monotone adversary.
        [ V, good_vector, tagged_vector ] = monotone_adversary( A, a, b );
        display(sum(tagged_vector))
        display(sum(sum(A-V)))
        display(sum(good_vector))
        % [ clustering_A, eigenvalues_A, eigenvectors_A ] = ...
        % spectral_clustering( A_store(:,:,repitition) );
        % [ clustering_V, eigenvalues_V, eigenvectors_V ] = ...
        % spectral_clustering( V_store(:,:,repitition) );
        [ clustering_A, eigenvalues_A, eigenvectors_A ] = spectral_clustering( A );
        [ clustering_V, eigenvalues_V, eigenvectors_V ] = spectral_clustering( V );
        error_rate_A = compute_error_rate(clustering_A,z_sbm);
        error_rate_V = compute_error_rate(clustering_V,z_sbm);
        disp(['Error rate for recovery from A is ' num2str(error_rate_A)])
        disp(['Error rate for recovery from V is ' num2str(error_rate_V)])
        if error_rate_A < (1 - percentage)
            disp(['GOOD RESULT on random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_A(iteration) = fail_records_A(iteration) + 1;
        end
        if error_rate_V < (1 - percentage)
            disp(['GOOD RESULT on semi-random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on semi-random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_V(iteration) = fail_records_V(iteration) + 1;
        end
        disp(['Finishing on repitition ' num2str(repitition) '.'])
    end
    disp(['Finishing on iteration ' num2str(iteration) '.'])
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Summary');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
for iteration = 1:num_of_trails
    jump = iteration/1;
    disp(['The trail failed ' num2str(fail_records_A(iteration)) ' times ' ...
    'on random model with gap ' num2str(jump) ' from threshold.']);
    disp(['The trail failed ' num2str(fail_records_V(iteration)) ' times ' ...
    'on semi-random model with gap ' num2str(jump) ' from threshold.']);
end
