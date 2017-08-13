% for i = 1:5
%     A = generate_sbm_adjacency(20,10,5);
%     V = monotone_adversary_random( A, 1/10 );
%     display(A)
%     display(V)
% end

rng('shuffle');
A = generate_sbm_adjacency( 1000, 5, 2 );
[ X, Xcost, D ] = sdp_solver(A);
% display(X)
% display(Xcost)
% display(D)
disp(['Check dual optimal: ' num2str(trace(D))]);
