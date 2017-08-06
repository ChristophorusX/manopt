# Autogenerated with SMOP 
from smop.core import *
# robust_pca.m

    
@function
def robust_pca(X=None,d=None,*args,**kwargs):
    varargin = robust_pca.varargin
    nargin = robust_pca.nargin

    # Computes a robust version of PCA (principal component analysis) on data.
# 
# function [U, cost] = robustpca(X, d)
    
    # Given a matrix X of size p by n, such that each column represents a
# point in R^p, this computes U: an orthonormal basis of size p by d such
# that the column space of U captures the points X as well as possible.
# More precisely, the function attempts to compute U as the minimizer
# over the Grassmann manifold (the set of linear subspaces) of:
    
    #  f(U) = (1/n) Sum_{i = 1:n} dist(X(:, i), the space spanned by U)
#       = (1/n) Sum_{i = 1:n} || U*U'*X(:, i) - X(:, i) ||
    
    # The output cost represents the average distance achieved with the
# returned U. Notice that norms are not squared, for robustness.
    
    # In practice, because this function is nonsmooth, it is smoothed with a
# pseudo-Huber loss function of parameter epsilon (noted e for short), and
# the smoothing parameter is iteratively reduced (with warm starts):
    
    #   f_e(U) = (1/n) Sum_{i = 1:n} l_e(|| U*U'*X(:, i) - X(:, i) ||)
    
    #   with l_e(x) = sqrt(x^2 + e^2) - e (for e = 0, this is absolute value).
    
    # The intermediate optimization of the smooth cost over the Grassmann
# manifold is performed using the Manopt toolbox.
    
    # Ideally, the non-outlier data should be centered. If not, this
# pre-processing centers all the data, but bear in mind that outliers will
# shift the center of mass too.
# X = X - repmat(mean(X, 2), [1, size(X, 2)]);
    
    # There are no guarantees that this code will return the optimal U.
# This code is distributed to illustrate one possible way of optimizing
# a nonsmooth cost function over a manifold, using Manopt with smoothing.
# For practical use, the constants in the code would need to be tuned.
    
    # This file is part of Manopt and is copyrighted. See the license file.
    
    # Main author: Nicolas Boumal and Teng Zhang, May 2, 2014
# Contributors:
    
    # Change log:
    
    #   March 4, 2015 (NB):
#       Uses a pseudo-Huber loss rather than a Huber loss: this has the
#       nice advantage of being smooth and simpler to code (no if's).
    
    #   April 8, 2015 (NB):
#       Built-in test data for quick tests; added comment about centering.
    
    # If no inputs, generate random data for illustration purposes.
    if nargin == 0:
        # Generate some data points aligned on a subspace
        X=dot(rand(2,1),(arange(1,30))) + multiply(dot(0.05,randn(2,30)),cat([(arange(1,30))],[(arange(1,30))]))
# robust_pca.m:58
        P=randperm(size(X,2))
# robust_pca.m:60
        outliers=10
# robust_pca.m:61
        X[:,P[1:outliers]]=dot(30,randn(2,outliers))
# robust_pca.m:62
        # X = X - repmat(mean(X, 2), [1, size(X, 2)]);
        d=1
# robust_pca.m:65
    
    # Prepare a Manopt problem structure for optimization of the given
    # cost (defined below) over the Grassmann manifold.
    p,n=size(X,nargout=2)
# robust_pca.m:74
    manifold=grassmannfactory(p,d)
# robust_pca.m:75
    problem.M = copy(manifold)
# robust_pca.m:76
    problem.cost = copy(robustpca_cost)
# robust_pca.m:77
    problem.egrad = copy(robustpca_gradient)
# robust_pca.m:78
    
    # This is just one idea: it is not necessarily useful or ideal.
    # Using a random initial guess, and starting over for a few different
    # ones is probably much better. For this example, we keep it simple.
    U,__,__=svds(X,d,nargout=3)
# robust_pca.m:84
    
    # the cost function over Grassmann.
    epsilon=1
# robust_pca.m:89
    n_iterations=6
# robust_pca.m:90
    reduction=0.5
# robust_pca.m:91
    options.verbosity = copy(2)
# robust_pca.m:92
    
    warning('off','manopt:getHessian:approx')
    for iter in arange(1,n_iterations).reshape(-1):
        U=trustregions(problem,U,options)
# robust_pca.m:95
        epsilon=dot(epsilon,reduction)
# robust_pca.m:96
    
    warning('on','manopt:getHessian:approx')
    
    epsilon=0
# robust_pca.m:102
    cost=robustpca_cost(U)
# robust_pca.m:103
    
    if nargin == 0:
        scatter(X[1,:],X[2,:])
        hold('on')
        plot(dot(dot(U[1],cat(- 1,1)),100),dot(dot(U[2],cat(- 1,1)),100),'r')
        hold('off')
        Upca,__,__=svds(X,1,nargout=3)
# robust_pca.m:114
        hold('on')
        plot(dot(dot(Upca[1],cat(- 1,1)),100),dot(dot(Upca[2],cat(- 1,1)),100),'k')
        hold('off')
        xlim(dot(1.1,cat(min(X[1,:]),max(X[1,:]))))
        ylim(dot(1.1,cat(min(X[2,:]),max(X[2,:]))))
        legend('data points','Robust PCA fit','Standard PCA fit')
    
    
    
    # Smoothed cost
    
@function
def robustpca_cost(U=None,*args,**kwargs):
    varargin = robustpca_cost.varargin
    nargin = robustpca_cost.nargin

    vecs=dot(U,(dot(U.T,X))) - X
# robust_pca.m:128
    sqnrms=sum(vecs ** 2,1)
# robust_pca.m:129
    vals=sqrt(sqnrms + epsilon ** 2) - epsilon
# robust_pca.m:130
    value=mean(vals)
# robust_pca.m:131
    return value
    
if __name__ == '__main__':
    pass
    
    # Euclidean gradient of the smoothed cost (it will be transformed into
    # the Riemannian gradient automatically by Manopt).
    
@function
def robustpca_gradient(U=None,*args,**kwargs):
    varargin = robustpca_gradient.varargin
    nargin = robustpca_gradient.nargin

    # Note that the computation of vecs and sqnrms is redundant
		# with their computation in the cost function. To speed
		# up the code, it would be wise to use the caching capabilities
		# of Manopt (the store structure). See online documentation.
		# It is not done here to keep the code a bit simpler.
    UtX=dot(U.T,X)
# robust_pca.m:144
    vecs=dot(U,UtX) - X
# robust_pca.m:145
    sqnrms=sum(vecs ** 2,1)
# robust_pca.m:146
    
    # and faster to compute the gradient.
        # G = zeros(p, d);
        # for i=1:n
        #     G = G + (1/sqrt(sqnrms(i) + epsilon^2)) * vecs(:,i) * UtX(:,i)';
        # end
        # G = G/n;
    G=mean(multiscale(1.0 / sqrt(sqnrms + epsilon ** 2),multiprod(reshape(vecs,cat(p,1,n)),multitransp(reshape(UtX,cat(d,1,n))))),3)
# robust_pca.m:154
    return G
    
if __name__ == '__main__':
    pass
    
    return G
    
if __name__ == '__main__':
    pass
    