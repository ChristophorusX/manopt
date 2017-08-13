function [ Y, z_syn ] = generate_synchronization_matrix( n, percent_of_elements_being_one, lambda )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random matrix Y of dimension n for Z2 synchronization problem,
% where Y=zz^T+\sigma W, with W a Wigner matrix generated from Gaussian
% distribution.
%
% @parameter: The number of elements in the group n.
% @parameter: The rough percentage percent_of_elements_being_one of vertices being in
% group 1.
% @parameter: The signal-to-noise ratio lambda.
%
% @return: The observed random matrix Y (with noise).
% @return: The ground truth of partition z_syn given in +1 and -1.
% NOTE: this setting is for synchronization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle');
z_pre = rand(n,1) > (1-percent_of_elements_being_one);
z_syn = (2*z_pre)-1;
% display(z_syn)
R = normrnd(0,1,n,n);
W = R - diag(R);
% display(W)
sigma = sqrt(n)/lambda; % lambda is the value to move around.
Y = z_syn*z_syn'+ sigma*W; % <= NOTE: Y is observed matrix of synchronization
% display(Y)
disp('Observed Z2 synchronization matrix is generated!')
end
