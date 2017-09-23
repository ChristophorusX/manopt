%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform GPM on the largest eigenvector of the observed matrix Y
% on Z2 synchronization.
%
% @author: Ruitu Xu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a synchronization matrix
n = 1000;
percentage = 0.5;
lambda = 10;
Y = generate_synchronization_gaussian(n, percentage, lambda);
% Generate linear part of the GPM
L = Y / n;
[eigen_vec, V] = eigs(Y);
largest = 1;
largest_val = 0;
for iter = 1:n
    if eigen_vec(iter) > largest_val
        disp('Encounter a larger eigenvalue!');
        largest = iter;
        disp(['Updated the largest one to be the ' str2num(iter) 'th eigenvalue.']);
        largest_val = eigen_vec(iter);
    end
end
x_0 = V(:,largest);
x = x_0;
for iter = 1:n
    l = L * x;
    x = 2 * (l > 0) - 1;
end
error_rate = compute_error_rate(x, z);
disp(['The error rate is ' error_rate]);
