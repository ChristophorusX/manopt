function [ X, Xcost, D_planted, eigenvalues, D_dual ] = sdp_solver_sym( Y, z )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a symmetric cost matrix Y and perform semidefinite
% programming under constraint X_{ii}=1 and X being positive semidefinite.
% @author: Ruitu Xu
%
% @parameter: Y be a symmetric cost function given by the problem setting.
% @parameter: z be the planted partion given by the generator of the model.
%
% @return: X be an optimum given by the sdp (symmetric).
% @return: Xcost be the optimal value of the objective function.
% @return: D_planted be the diagonal matrix used in spectral analysis.
% @return: eigenvalues be the list of eigenvalues needed in spectral
%          analysis (ordered in ascending order).
% @return: D_dual be dual variable produced by sdp.
%
% Note: this function utilizes cvx package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ dimension, ~ ] = size(Y);

cvx_clear

cvx_begin sdp
variable X(dimension,dimension) symmetric;
dual variable D;
maximize( trace(Y*X) );
subject to
    D : diag(X) == 1;
    X >= 0;
cvx_end

Xcost = trace(Y*X);
D_dual = -diag(D);

D_planted_pre = diag(z) * Y * diag(z);
D_planted = diag(sum(D_planted_pre, 2));
e = eig(D_planted - B);
eigenvalues_pre = sort(e);
eigenvalues = eigenvalues_pre .* (eigenvalues_pre > eps)

end  % function
