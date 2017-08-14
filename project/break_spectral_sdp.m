%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters and trying to break the spectral method
% and sdp method on it (in logrithmic regime).
% @author: Ruitu Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1000;
a = 25;
num_of_trails = 10;
num_of_repititions = 1;
density_of_jump = 1;
delta = 1;
b_base = 10;
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
        [ clustering_A, eigenvalues_A, eigenvectors_A ] = spectral_clustering( A );
        [ clustering_V_MPW, eigenvalues_V_MPW, eigenvectors_V_MPW ] = spectral_clustering( V_MPW );
        [ clustering_V_hub, eigenvalues_V_hub, eigenvectors_V_hub ] = spectral_clustering( V_hub );
        [ clustering_V_random, eigenvalues_V_random, eigenvectors_V_random ] = spectral_clustering( V_random );
        error_rate_A = compute_error_rate(clustering_A,z_sbm);
        error_rate_V_MPW = compute_error_rate(clustering_V_MPW,z_sbm);
        error_rate_V_hub = compute_error_rate(clustering_V_hub,z_sbm);
        error_rate_V_random = compute_error_rate(clustering_V_random,z_sbm);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Start printing evaluation output on Spectral-SBM.')
        disp('=========================================================================')
        disp(['Error rate for recovery from A is ' num2str(error_rate_A)])
        disp(['Error rate for recovery from V_MPW is ' num2str(error_rate_V_MPW)])
        disp(['Error rate for recovery from V_hub is ' num2str(error_rate_V_hub)])
        disp(['Error rate for recovery from V_random is ' num2str(error_rate_V_random)])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Exact recovery??????')
        if error_rate_A == 0
            disp('Yes, it achieves exact recovery on random model.');
        else
            disp('No, it does not achieve exact recovery on random model.');
        end
        if error_rate_V_MPW == 0
            disp('Yes, it achieves exact recovery on semi-random model under MPW adversary.');
        else
            disp('No, it does not achieve exact recovery on semi-random model under MPW adversary.');
        end
        if error_rate_V_hub == 0
            disp('Yes, it achieves exact recovery on semi-random model under hub adversary.');
        else
            disp('No, it does not achieve exact recovery on semi-random model under hub adversary.');
        end
        if error_rate_V_random == 0
            disp('Yes, it achieves exact recovery on semi-random model under random adversary.');
        else
            disp('No, it does not achieve exact recovery on semi-random model under random adversary.');
        end
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        all_ones = ones(n, 1);
        B = 2 * A - (all_ones * all_ones' + eye(n));
        B_MPW = 2 * V_MPW - (all_ones * all_ones' + eye(n));
        B_hub = 2 * V_hub - (all_ones * all_ones' + eye(n));
        B_random = 2 * V_random - (all_ones * all_ones' + eye(n));
        [ ~, Xcost_B, ~, eigenvalues_B, D_dual_B ] = sdp_solver_sym(B, z_sbm);
        [ ~, Xcost_B_MPW, ~, eigenvalues_B_MPW, D_dual_B_MPW ] = sdp_solver_sym(B_MPW, z_sbm);
        [ ~, Xcost_B_hub, ~, eigenvalues_B_hub, D_dual_B_hub ] = sdp_solver_sym(B_hub, z_sbm);
        [ ~, Xcost_B_random, ~, eigenvalues_B_random, D_dual_B_random ] = sdp_solver_sym(B_random, z_sbm);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Start printing evaluation output on SDP-SBM.')
        disp('=========================================================================')
        disp('The smallest 3 eigenvalues are: ');
        for subsubiter = 1:3
            disp(num2str(eigenvalues(subsubiter)));
        end
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Exact recovery??????')
        if eigenvalues_B(2) > eps
            disp('Yes, it achieves exact recovery on random model.');
        else
            disp('No, it does not achieve exact recovery on random model.');
        end
        if eigenvalues_B_MPW(2) > eps
            disp('Yes, it achieves exact recovery on semi-random model under MPW adversary.');
        else
            disp('No, it does not achieve exact recovery on semi-random model under MPW adversary.');
        end
        if eigenvalues_B_hub(2) > eps
            disp('Yes, it achieves exact recovery on semi-random model under MPW adversary.');
        else
            disp('No, it does not achieve exact recovery on semi-random model under MPW adversary.');
        end
        if eigenvalues_B_random(2) > eps
            disp('Yes, it achieves exact recovery on semi-random model under MPW adversary.');
        else
            disp('No, it does not achieve exact recovery on semi-random model under MPW adversary.');
        end
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
end
