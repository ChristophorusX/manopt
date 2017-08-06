% This is a MATLAB program of running monotone adversary on community
% detection problem of 2 clusters. @author: Ruitu Xu

% Generate random adjacency matrix of dimension n (even) according to inner and
% inter cluster edge probability p and q
number_of_vertices = 1000;
a = 4;
b = 6;
n = number_of_vertices;
p = a * log(n) / n;
q = b * log(n) / n;
R1 = rand(n/2,n/2) > (1-p);
R2 = rand(n/2,n/2) > (1-q);
A = [R1 R2; R2 R1];
row_sum = sum(A,2);
for i = 1:n
    disp(['The summation of row ' num2str(i) ' is ' num2str(row_sum(i)) '.'])
end
disp(['The average edge of every row should be around ' num2str((p*n + q*n)/2) '.'])
disp(['The empirical average number of edges is ' num2str(sum(row_sum)/n)])
    

