%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a program verifying the spectral gap increase of synchronization
% problem under rank 2 Burer-Monteiro approach against monotone adversary.
%
% @author: Ruitu Xu
% NOTE: this setting is for synchronization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'manopt:getHessian:approx')
% Generating the observed advhronization matrix.
% Generating the advhronization matrix under monotone adversary.
[ Y_sync, Y_adv ] = monotone_adversary_sync(100, 0.5, 10, 1);
% disp(((Y_adv - Y_sync) .* Y_sync) > 0);
% Perform BM on both matrix.
[ Q_Y_sync, ~, ~, ~ ] = burer_monteiro( Y_sync );
[ Q_Y_adv, ~, ~, ~ ] = burer_monteiro( Y_adv );
% Construct S
S_sync = diag(diag(Y_sync * (Q_Y_sync * Q_Y_sync'))) - Y_sync;
S_adv = diag(diag(Y_adv * (Q_Y_adv * Q_Y_adv'))) - Y_adv;
% Compute eigenvalues of S
e_sync = eig(S_sync);
e_adv = eig(S_adv);
% Sorting eigenvalues
eigenvalues_sync = sort(e_sync);
eigenvalues_adv = sort(e_adv);
output_sync = 'spectrum of the sync is: ';
output_adv = 'spectrum of the adv is: ';
output_spectral_gap_increase = 'spectral gap increase is: ';
for iter = 1:100
    output_sync = strcat(output_sync, num2str(eigenvalues_sync(iter)), ', ');
    output_adv = strcat(output_adv, num2str(eigenvalues_adv(iter)), ', ');
    if eigenvalues_adv(iter) >= eigenvalues_sync(iter)
        output_spectral_gap_increase = strcat(output_spectral_gap_increase, 'LARGER', ', ');
    else
        output_spectral_gap_increase = strcat(output_spectral_gap_increase, 'smaller', ', ');
    end
end
disp(output_sync);
disp(output_adv);
% Compute spectral gap increase.
disp(output_spectral_gap_increase);
