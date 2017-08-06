# Autogenerated with SMOP 
from smop.core import *
# maxcut.m

    
@function
def maxcut(L=None,r=None,*args,**kwargs):
    varargin = maxcut.varargin
    nargin = maxcut.nargin

    # Algorithm to (try to) compute a maximum cut of a graph, via SDP approach.
# 
# function x = maxcut(L)
# function [x, cutvalue, cutvalue_upperbound, Y] = maxcut(L, r)
    
    # L is the Laplacian matrix describing the graph to cut. The Laplacian of a
# graph is L = D - A, where D is the diagonal degree matrix (D(i, i) is the
# sum of the weights of the edges adjacent to node i) and A is the
# symmetric adjacency matrix of the graph (A(i, j) = A(j, i) is the weight
# of the edge joining nodes i and j). If L is sparse, this will be
# exploited.
    
    # If the graph has n nodes, then L is nxn and the output x is a vector of
# length n such that x(i) is +1 or -1. This partitions the nodes of the
# graph in two classes, in an attempt to maximize the sum of the weights of
# the edges that go from one class to the other (MAX CUT problem).
    
    # cutvalue is the sum of the weights of the edges 'cut' by the partition x.
    
    # If the algorithm reached the global optimum of the underlying SDP
# problem, then it produces an upperbound on the maximum cut value. This
# value is returned in cutvalue_upperbound if it is found. Otherwise, that
# output is set to NaN.
    
    # If r is specified (by default, r = n), the algorithm will stop at rank r.
# This may prevent the algorithm from reaching a globally optimal solution
# for the underlying SDP problem (but can greatly help in keeping the
# execution time under control). If a global optimum of the SDP is reached
# before rank r, the algorithm will stop of course.
    
    # Y is a matrix of size nxp, with p <= r, such that X = Y*Y' is the best
# solution found for the underlying SDP problem. If cutvalue_upperbound is
# not NaN, then Y*Y' is optimal for the SDP and cutvalue_upperbound is its
# cut value.
# 
# By Goemans and Williamson 1995, it is known that if the optimal value of
# the SDP is reached, then the returned cut, in expectation, is at most at
# a fraction 0.878 of the optimal cut. (This is not exactly valid because
# we do not use random projection here; sign(Y*randn(size(Y, 2), 1)) will
# give a cut that respects this statement -- it's usually worse though).
    
    # The algorithm is essentially that of:
# Journee, Bach, Absil and Sepulchre, SIAM 2010
# Low-rank optimization on the cone of positive semidefinite matrices.
    
    # It is itself based on the famous SDP relaxation of MAX CUT:
# Goemans and Williamson, 1995
# Improved approximation algorithms for maximum cut and satisfiability
# problems using semidefinite programming.
# 
# See also: elliptope_SDP
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, July 18, 2013
# Contributors:
# Change log:
#   
#   April 3, 2015 (NB):
#       L products now counted with the new shared memory system. This is
#       more reliable and more flexible than using a global variable.
    
    # If no inputs are provided, generate a random graph Laplacian.
    # This is for illustration purposes only.
    if logical_not(exist('L','var')) or isempty(L):
        n=20
# maxcut.m:67
        A=triu(randn(n) <= 0.4,1)
# maxcut.m:68
        A=A + A.T
# maxcut.m:69
        D=diag(sum(A,2))
# maxcut.m:70
        L=D - A
# maxcut.m:71
    
    n=size(L,1)
# maxcut.m:75
    assert_(size(L,2) == n,'L must be square.')
    if logical_not(exist('r','var')) or isempty(r) or r > n:
        r=copy(n)
# maxcut.m:79
    
    
    # We will let the rank increase. Each rank value will generate a cut.
    # We have to go up in the rank to eventually find a certificate of SDP
    # optimality. This in turn will provide an upperbound on the MAX CUT
    # value and ensure that we're doing well, according to Goemans and
    # Williamson's argument. In practice though, the good cuts often come
    # up for low rank values, so we better keep track of the best one.
    best_x=ones(n,1)
# maxcut.m:88
    best_cutvalue=0
# maxcut.m:89
    cutvalue_upperbound=copy(NaN)
# maxcut.m:90
    time=matlabarray([])
# maxcut.m:92
    cost=matlabarray([])
# maxcut.m:93
    for rr in arange(2,r).reshape(-1):
        manifold=elliptopefactory(n,rr)
# maxcut.m:97
        if rr == 2:
            # At first, for rank 2, generate a random point.
            Y0=manifold.rand()
# maxcut.m:102
        else:
            # To increase the rank, we could just add a column of zeros to
            # the Y matrix. Unfortunately, this lands us in a saddle point.
            # To escape from the saddle, we may compute an eigenvector of
            # Sy associated to a negative eigenvalue: that will yield a
            # (second order) descent direction Z. See Journee et al ; Sy is
            # linked to dual certificates for the SDP.
            Y0=matlabarray(cat(Y,zeros(n,1)))
# maxcut.m:112
            LY0=dot(L,Y0)
# maxcut.m:113
            Dy=spdiags(sum(multiply(LY0,Y0),2),0,n,n)
# maxcut.m:114
            Sy=(Dy - L) / 4
# maxcut.m:115
            eigsopts.issym = copy(true)
# maxcut.m:117
            eigsopts.isreal = copy(true)
# maxcut.m:118
            v,s=eigs(Sy,1,'SA',eigsopts,nargout=2)
# maxcut.m:119
            # a saddle point: we're actually done!
            if s >= - 1e-08:
                # We can stop here: we found the global optimum of the SDP,
                # and hence the reached cost is a valid upper bound on the
                # maximum cut value.
                cutvalue_upperbound=max(- cat(info.cost))
# maxcut.m:126
                break
            # This is our escape direction.
            Z=manifold.proj(Y0,cat(zeros(n,rr - 1),v))
# maxcut.m:131
            # # function looks like at a saddle point. But will require the
            # # problem structure which is not defined here: see the helper
            # # function.
            # plotprofile(problem, Y0, Z, linspace(-1, 1, 101));
            # drawnow; pause;
            # Now make a step in the Z direction to escape from the saddle.
            # It is not obvious that it is ok to do a unit step ... perhaps
            # need to be cautious here with the stepsize. It's not too
            # critical though: the important point is to leave the saddle
            # point. But it's nice to guarantee monotone decrease of the
            # cost, and we can't do that with a constant step (at least,
            # not without a proper argument to back it up).
            stepsize=1
# maxcut.m:147
            Y0=manifold.retr(Y0,Z,stepsize)
# maxcut.m:148
        # Use the Riemannian optimization based algorithm lower in this
        # file to reach a critical point (typically a local optimizer) of
        # the max cut cost with fixed rank, starting from Y0.
        Y,info=maxcut_fixedrank(L,Y0,nargout=2)
# maxcut.m:155
        thistime=matlabarray(cat(info.time))
# maxcut.m:158
        if logical_not(isempty(time)):
            thistime=time[end()] + thistime
# maxcut.m:160
        time=matlabarray(cat(time,thistime))
# maxcut.m:162
        cost=matlabarray(cat(cost,cat(info.cost)))
# maxcut.m:163
        # Time to turn the matrix Y into a cut.
        # We can either do the random rounding as follows:
        # x = sign(Y*randn(rr, 1));
        # or extract the "PCA direction" of the points in Y and cut
        # orthogonally to that direction, as follows (seems faster than
        # calling svds):
        U,__,__=svd(Y,0,nargout=3)
# maxcut.m:171
        u=U[:,1]
# maxcut.m:172
        x=sign(u)
# maxcut.m:173
        cutvalue=(dot(dot(x.T,L),x)) / 4
# maxcut.m:175
        if cutvalue > best_cutvalue:
            best_x=copy(x)
# maxcut.m:177
            best_cutvalue=copy(cutvalue)
# maxcut.m:178
    
    
    x=copy(best_x)
# maxcut.m:183
    cutvalue=copy(best_cutvalue)
# maxcut.m:184
    plot(time,- cost,'.-')
    xlabel('Time [s]')
    ylabel('Relaxed cut value')
    title('The relaxed cut value is an upper bound on the optimal cut value.')
    return x,cutvalue,cutvalue_upperbound,Y
    
if __name__ == '__main__':
    pass
    
    
@function
def maxcut_fixedrank(L=None,Y=None,*args,**kwargs):
    varargin = maxcut_fixedrank.varargin
    nargin = maxcut_fixedrank.nargin

    # Try to solve the (fixed) rank r relaxed max cut program, based on the
# Laplacian of the graph L and an initial guess Y. L is nxn and Y is nxr.
    
    n,r=size(Y,nargout=2)
# maxcut.m:198
    assert_(all(size(L) == n))
    
    # semidefinite matrices of size n with rank r and all diagonal entries
    # are 1.
    manifold=elliptopefactory(n,r)
# maxcut.m:204
    
    # # against the (conceptually simpler) oblique manifold geometry,
    # # uncomment this line.
    # manifold = obliquefactory(r, n, true);
    
    problem.M = copy(manifold)
# maxcut.m:211
    
    # # function and its gradient and Hessian (here expressed using the
    # # Euclidean gradient and Hessian).
    # problem.cost  = @(Y)  -trace(Y'*L*Y)/4;
    # problem.egrad = @(Y) -(L*Y)/2;
    # problem.ehess = @(Y, U) -(L*U)/2;
    
    # Instead of the prototyping version, the functions below describe the
    # cost, gradient and Hessian using the caching system (the store
    # structure). This alows to execute exactly the required number of
    # multiplications with the matrix L. These multiplications are counted
    # using the shared memory in the store structure: that memory is
    # shared , so we get access to the same data, regardless of the
    # point Y currently visited.
    
    # For every visited point Y, we will need L*Y. This function makes sure
    # the quantity L*Y is available, but only computes it if it wasn't
    # already computed.
    
@function
def prepare(Y=None,store=None,*args,**kwargs):
    varargin = prepare.varargin
    nargin = prepare.nargin

    if logical_not(isfield(store,'LY')):
        # Compute and store the product for the current point Y.
        store.LY = copy(dot(L,Y))
# maxcut.m:234
        if isfield(store.shared,'counter'):
            store.shared.counter = copy(store.shared.counter + 1)
# maxcut.m:237
        else:
            store.shared.counter = copy(1)
# maxcut.m:239
    
    return store
    
if __name__ == '__main__':
    pass
    
    problem.cost = copy(cost)
# maxcut.m:244
    
@function
def cost(Y=None,store=None,*args,**kwargs):
    varargin = cost.varargin
    nargin = cost.nargin

    store=prepare(Y,store)
# maxcut.m:246
    LY=store.LY
# maxcut.m:247
    f=- (dot(ravel(Y).T,ravel(LY))) / 4
# maxcut.m:248
    
    return f,store
    
if __name__ == '__main__':
    pass
    
    problem.egrad = copy(egrad)
# maxcut.m:251
    
@function
def egrad(Y=None,store=None,*args,**kwargs):
    varargin = egrad.varargin
    nargin = egrad.nargin

    store=prepare(Y,store)
# maxcut.m:253
    LY=store.LY
# maxcut.m:254
    g=- LY / 2
# maxcut.m:255
    return g,store
    
if __name__ == '__main__':
    pass
    
    problem.ehess = copy(ehess)
# maxcut.m:258
    
@function
def ehess(Y=None,U=None,store=None,*args,**kwargs):
    varargin = ehess.varargin
    nargin = ehess.nargin

    store=prepare(Y,store)
# maxcut.m:260
    
    LU=dot(L,U)
# maxcut.m:261
    store.shared.counter = copy(store.shared.counter + 1)
# maxcut.m:262
    h=- LU / 2
# maxcut.m:263
    return h,store
    
if __name__ == '__main__':
    pass
    
    # statsfun is called exactly once after each iteration (including after
    # the evaluation of the cost at the initial guess). We then register
    # the value of the L-products counter (which counts how many products
    # were needed so far).
    # options.statsfun = @statsfun;
    # function stats = statsfun(problem, Y, stats, store) ##ok
    #     stats.Lproducts = store.shared.counter;
    # end
    # Equivalent, but simpler syntax:
    options.statsfun = copy(statsfunhelper('Lproducts',lambda problem=None,Y=None,stats=None,store=None: store.shared.counter))
# maxcut.m:275
    
    # # correct during the prototyping stage.
    # checkgradient(problem); pause;
    # checkhessian(problem); pause;
    
    # # To investigate the effect of the rotational invariance when using
    # # the oblique or the elliptope geometry, or to study the saddle point
    # # issue mentioned above, it is sometimes interesting to look at the
    # # spectrum of the Hessian. For large dimensions, this is slow!
    # stairs(sort(hessianspectrum(problem, Y)));
    # drawnow; pause;
    
    
    # # When facing a saddle point issue as described in the master
    # # function, and when no sure mechanism exists to find an escape
    # # direction, it may be helpful to set useRand to true and raise
    # # miniter to more than 1, when using trustregions. This will tell the
    # # solver to not stop before at least miniter iterations were
    # # accomplished (thus disregarding the zero gradient at the saddle
    # # point) and to use random search directions to kick start the inner
    # # solve (tCG) step. It is not as efficient as finding a sure escape
    # # direction, but sometimes it's the best we have.
    # options.useRand = true;
    # options.miniter = 5;
    
    options.verbosity = copy(2)
# maxcut.m:304
    Y,Ycost,info=trustregions(problem,Y,options,nargout=3)
# maxcut.m:305
    
    
    fprintf('Products with L: %d\\n',max(cat(info.Lproducts)))
    return h,store
    
if __name__ == '__main__':
    pass
    