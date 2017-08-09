function [ Q, Qcost, info, options ] = burer_monteiro( B )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform a manifold optimization on the cost function Tr(Q^TBQ) given the
% manifold (Symmetric positive semidefinite, fixed-rank with unit diagonal)
% with rank(X=QQ^T)<=2.
% It receives
% - A matrix to perform Burer-Monteiro on.
% It returns
% - The resulting n by 2 matrix Q.
% - The value of cost function under Q.
% - Information returned by the trust region method.
% - Options used by the trust region method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2;
dimensions = size(B);
dimension = dimensions(1);
manifold = elliptopefactory(dimension, k);
problem.M = manifold;
% Note that Manopt only minimize the cost function.
problem.cost = @(Q) -trace(Q'*B*Q);
% The Euclidian gradient is given by calculation.
problem.egrad = @(Q) -((B+B')*Q+[diag(B*diag(Q(:,1))) diag(B*diag(Q(:,2)))]);
% Numerically check gradient and Hessian consistency.
% figure;
% checkgradient(problem);
% figure;
% checkhessian(problem);
% Pass the problem to the solver, i.e. trust regions.
[Q, Qcost, info, options] = trustregions(problem);
% disp(['The output optimal cost by Burer-Monteiro is ' num2str(Qcost) '.'])

% Figure on convergence of the algorithm.
% figure;
% semilogy([info.iter], [info.gradnorm], '.-');
% xlabel('Iteration #');
% ylabel('Gradient norm');
% title('Convergence of the trust-regions algorithm on the manifold');


end

