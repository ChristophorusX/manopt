function [ Q, Qcost, info, options ] = burer_monteiro_k( B )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform a manifold optimization on the cost function Tr(Q^TBQ) given the
% manifold (Symmetric positive semidefinite, fixed-rank with unit diagonal)
% with rank(X=QQ^T)<=2.
% @author: Ruitu Xu
%
% @parameter: A matrix B to perform Burer-Monteiro on.
%
% @return: The resulting n by 2 matrix Q.
% @return: The value of cost function under Q.
% @return: Information returned by the trust region method.
% @return: Options used by the trust region method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ dimension, ~ ] = size(B);
k = round(sqrt(2*dimension));
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
disp('Burer-Monteiro approach finished!')

end
