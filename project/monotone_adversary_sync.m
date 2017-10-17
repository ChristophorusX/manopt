function [ Y_sync, Y_adv ] = monotone_adversary_sync( n, percent_of_elements_being_one, lambda, sigma )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random matrix Y_adv of dimension n for Z2 synchronization
% problem by modifying the observed matrix with a Gaussian monotone
% noise.
% @author: Ruitu Xu
%
% @parameter: Dimension n.
% @parameter: percent_of_elements_being_one.
% @parameter: SNR lambda.
% @parameter: The standard deviation sigma of the monotone noise.
%
% @return: The observed random matrix Y_sync.
% @return: The observed random matrix Y_adv (with monotone noise).
% NOTE: this setting is for synchronization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle');
z_pre = rand(n,1) > (1-percent_of_elements_being_one);
z_syn = (2*z_pre)-1;
% display(z_syn)
R_pre = normrnd(0,sqrt(2)/2,n,n);
R = R_pre + R_pre'; % Sum of independent normal is also normal
W = R - diag(R);
% display(W)
sigma = sqrt(n)/lambda; % lambda is the value to move around.
Y_sync = z_syn*z_syn'+ sigma*W; % <= NOTE: Y is observed matrix of synchronization
disp('Observed Z2 synchronization matrix is generated!')

M = normrnd(0,sigma,n,n);
M_plus = abs(M);
M_pre = M_plus .* (((W > 0) + (W < 0)) > 0);
multiplier = (W > 0) * 2 - ones(n,n);
Y_adv = Y_sync + multiplier .* M_pre;
disp('Monotone adversary matrix of Z2 synchronization is generated!')
end
