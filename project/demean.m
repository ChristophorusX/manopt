function [ A_bm ] = demean( A, p, q )
% function: Short description
%
% Extended description

[ dimension, ~ ] = size( A );
one_vector = ones(dimension, 1);
A_bm = A - (p + q)/2 * (one_vector * one_vector');

end  % function
