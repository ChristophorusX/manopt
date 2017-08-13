function [ X, Xcost, D, eigenvalues ] = sdp_solver_sym( Y )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a symmetric cost matrix Y and perform semidefinite
% programming under constraint X_{ii}=1 and X being positive semidefinite.
% @author: Ruitu Xu
%
% @parameter: Y be a symmetric cost function given by the problem setting.
%
% @return: X be an optimum given by the sdp (symmetric).
% @return: Xcost be the optimal value of the objective function.
% @return: D be dual variable produced by sdp.
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

D = -diag(D);
Xcost = trace(Y*X);
e = eig(D - Y);
eigenvalues = sort(e);

end  % function
