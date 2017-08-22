function [ V_bm ] = demean_adversary( V )
% function: Short description
%
% Extended description

[ dimension, ~ ] = size(V);
sum_col = sum(V, 2);
one_vector = ones(dimension, 1);
V_bm = V - sum_col * one_vector';

end  % function
