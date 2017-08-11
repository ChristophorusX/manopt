%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and trying to break the spectral method
% on it (in sparse regime).
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix a percentage of partial recovery.
percentage = 0.6;
% Set a to a fixed number and the set b according to the treshold.
a = 20;
% Set an array recording the times of failing trails.
num_of_trails = 100;
fail_records_A = zeros(num_of_trails,1);
fail_records_V = zeros(num_of_trails,1);
for iteration = 1:num_of_trails
    jump = iteration/10;
    b = vpasolve((a-x)^2==(2+jump)*(a+x),x); % obtain b from a and threshold.
    % Generate adjacency matrix.
    for i=1:100
        [ A, z_sbm ] = generate_sbm_adjacency(n, a, b);
        % Put matrix A into monotone adversary.
        [ V, good_vector, tagged_vector ] = monotone_adversary( A, a, b );
        % display(sum(tagged_vector))
        % display(sum(sum(A-V)))
        % display(good_vector)
        [ clustering_A, eigenvalues_A, eigenvectors_A ] = spectral_clustering( A );
        [ clustering_V, eigenvalues_V, eigenvectors_V ] = spectral_clustering( V );
        error_rate_A = compute_error_rate(clustering_A,z_sbm);
        error_rate_V = compute_error_rate(clustering_V,z_sbm);
        if error_rate_A < (1 - percentage)
            disp(['GOOD RESULT on random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_A(iteration) += 1;
        end
        if error_rate_V < (1 - percentage)
            disp(['GOOD RESULT on semi-random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
        else
            disp(['BAD RESULT on semi-random model with parameter a = ' num2str(a) ', b = ' ...
            num2str(b) ', and over threshold with amount ' num2str(jump)]);
            fail_records_V(iteration) += 1;
        end
    end
end
