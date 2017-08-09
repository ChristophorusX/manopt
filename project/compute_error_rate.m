function [ error_rate ] = compute_error_rate( clustering, z )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the error rate given a clustering result and the ground truth.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimensions = size(z);
dimension = dimensions(1);
diff_1 = z - clustering;
diff_2 = z + clustering;
errors_1 = sum(diff_1 ~=0);
errors_2 = sum(diff_2 ~=0);
error_rate = min(errors_1,errors_2)/dimension;

end

