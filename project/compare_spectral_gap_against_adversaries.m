%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a program verifying the spectral gap increase of advhronization
% problem under rank 2 Burer-Monteiro approach against monotone adversary.
%
% @author: Ruitu Xu
% NOTE: this setting is for advhronization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'manopt:getHessian:approx')
% Generating the observed advhronization matrix.
Y_sync = generate_synchronization_gaussian(100, 0.5, 10);
disp(Y_sync);
% Generating the advhronization matrix under monotone adversary.
Y_adv = monotone_adversary_sync(Y_sync, 1);
disp(Y_adv);
% Perform BM on both matrix.
[ Q_Y_sync, ~, ~, ~ ] = burer_monteiro( Y_sync );
[ Q_Y_adv, ~, ~, ~ ] = burer_monteiro( Y_adv );
% Construct S
S_sync = diag(diag(Y_sync * Q_Y_sync * Q_Y_sync')) - Y_sync;
S_adv = diag(diag(Y_adv * Q_Y_adv * Q_Y_adv')) - Y_adv;
% Compute eigenvalues of S
e_sync = eig(S_sync);
e_adv = eig(S_adv);
% Sorting eigenvalues
eigenvalues_sync = sort(e_sync);
eigenvalues_adv = sort(e_adv);
% Compute spectral gap increase.
spectral_gap_increase = eigenvalues_adv(2) - eigenvalues_sync(2);
disp(['The spectral gap increase is ' num2str(spectral_gap_increase)]);
