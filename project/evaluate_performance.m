function [ true_cost_value, correlation ] = evaluate_performance( z, B, Q )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing some evaluations on the Burer-Monteiro methods on respective
% problems.
% @author: Ruitu Xu
%
% @parameter: A ground truth vector z.
% @parameter: A cost matrix of the problem B.
% @parameter: The second order critical point Q from Burer-Monteiro.
%
% @return: The cost value of the ground truth vector on cost matrix.
% @return: The correlation between ground truth z and critical point Q.
% The function also plots the row vectors of second order critical point Q
% on a unit circle to reveal their relative position and potential
% groupings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting evaluating the performance...')
% disp(z)
true_cost_value = -z'*B*z; % Cost given ground truth z.
% disp(true_cost_value)
% disp(['The planted cost value of this problem is ' num2str(true_cost_value) '.'])
% Compute the correlation between the critical point Q and ground truth z.
dimensions = size(B);
dimension = dimensions(1);
correlation = sqrt(z'*(Q*Q')*z)/dimension;
% disp(['The correlation between output X=QQ^T and the planted vector z is ' ...
%     num2str(correlation) '.'])

% figure;
% plot_x = Q(:,1);
% display(plot_x)
% plot_y = Q(:,2);
% display(plot_y)
% title('Row vectors of second order critical point Q on unit circle.');
% scatter(plot_x,plot_y)

end
