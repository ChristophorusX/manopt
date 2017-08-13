function [ X, Xcost, D ] = sdp_solver( Y )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a cost matrix Y and perform semidefinite programming
% under constraint X_{ii}=1 and X being positive semidefinite.
%
% @parameter: Y be a cost function given by the problem setting.
%
% @return: X be an optimum given by the sdp.
% @return: Xcost be the optimal value of the objective function.
% @return: D be dual variable produced by sdp.
%
% Note: this function utilizes cvx package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ dimension, ~ ] = size(Y);

cvx_clear

cvx_begin sdp
variable X(dimension,dimension) semidefinite;
dual variable D;
maximize( trace(Y*X) );
subject to
    D : diag(X) == 1;
    X >= 0;
cvx_end

D = diag(D);
Xcost = trace(Y*X);

end  % function
