% for i = 1:5
%     A = generate_sbm_adjacency(20,10,5);
%     V = monotone_adversary_random( A, 1/10 );
%     display(A)
%     display(V)
% end

rng('shuffle');
A = generate_sbm_adjacency_logrithmic( 1000, 5, 2 );
[ X, Xcost, D, eigenvalues ] = sdp_solver_sym(A);
disp(['Check dual optimal: ' num2str(trace(D))]);
disp('The smallest 10 eigenvalues are: ');
for iter = 1:10
    disp(num2str(eigenvalues(iter)));
end

[ X, Xcost, D, eigenvalues ] = sdp_solver(A);
disp(['Check dual optimal: ' num2str(trace(D))]);
disp('The smallest 10 eigenvalues are: ');
for iter = 1:10
    disp(num2str(eigenvalues(iter)));
end
% display(X)
% display(Xcost)
% display(D)
