function [ Y_adv ] = monotone_adversary_sync( Y_sync, sigma )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random matrix Y_adv of dimension n for Z2 synchronization
% problem by modifying the observed matrix with a Gaussian monotone
% noise.
% @author: Ruitu Xu
%
% @parameter: Observed synchronization matrix Y_sync.
% @parameter: The standard deviation sigma of the monotone noise.
%
% @return: The observed random matrix Y_sync (with monotone noise).
% NOTE: this setting is for synchronization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle');
[n,~] = size(Y_sync);
M = normrnd(0,sigma,n,n);
M_plus = abs(M);
M_pre = M_plus .* (((Y_sync > 0) + (Y_sync < 0)) > 0);
multiplier = (Y_sync > 0) * 2 - ones(n,n);
Y_adv = Y_sync + multiplier .* M_pre;
disp('Monotone adversary matrix of Z2 synchronization is generated!')
end
