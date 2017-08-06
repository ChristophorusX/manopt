# Autogenerated with SMOP 
from smop.core import *
# hessianextreme.m

    
@function
def hessianextreme(problem=None,x=None,side=None,y0=None,options=None,storedb=None,key=None,*args,**kwargs):
    varargin = hessianextreme.varargin
    nargin = hessianextreme.nargin

    # Compute an extreme eigenvector / eigenvalue of the Hessian of a problem.
    
    # [u, lambda, info] = hessianextreme(problem, x)
# [u, lambda, info] = hessianextreme(problem, x, side)
# [u, lambda, info] = hessianextreme(problem, x, side, u0)
# [u, lambda, info] = hessianextreme(problem, x, side, u0, options)
# [u, lambda, info] = hessianextreme(problem, x, side, u0, options, storedb)
# [u, lambda, info] = hessianextreme(problem, x, side, u0, options, storedb, key)
# 
# (For side, u0 and options, pass [] to omit any.)
    
    # Given a Manopt problem structure and a point x on the manifold problem.M,
# this function computes a tangent vector u at x of unit norm such that the
# Hessian quadratic form is minimized or maximized:
    
    #    minimize or maximize <u, Hess f(x)[u]> such that <u, u> = 1,
    
    # where <.,.> is the Riemannian metric on the tangent space at x. Choose
# between minimizing and maximizing by setting side = 'min' or 'max', with
# 'min' being the default. The value attained is returned as lambda, and
# is the minimal or maximal eigenvalue of the Hessian (actually, the last
# value attained when the solver stopped). This is a real number since the
# Hessian is a symmetric operator.
    
    # If u0 is specified, it should be a unit-norm tangent vector at x. It is
# then used as initial guess to solve the above problem. Pass [] to omit.
    
    # The options structure, if provided, will be passed along to manoptsolve.
# As such, you may choose which solver to use to solve the above
# optimization problem by setting options.solver. See manoptsolve's help.
# The other options will be passed along to the chosen solver too.
# Pass [] to omit.
    
    # Often times, it is only necessary to compute a vector u such that the
# quadratic form is negative, if that is at all possible. To do so, set the
# following stopping criterion: options.tolcost = -1e-10; (for example)
# and side = 'min'. The solver will return as soon as the quadratic cost
# defined above drops below the set value (or sooner if another stopping
# criterion triggers first.)
    
    # storedb is a StoreDB object, key is the StoreDB key to point x.
    
    # info is the info struct-array returned by the solver.
    
    # See also: hessianspectrum manoptsolve tangentspherefactory
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Aug. 13, 2014.
# Contributors: 
# Change log:
    
    #   April 3, 2015 (NB):
#       Works with the new StoreDB class system.
    
    #   May 7, 2015 (NB):
#       Default solver options: verbosity = 0 and defaults to trustregions.
    
    #   Nov 27, 2015 (NB):
#       The function now also returns the info struct-array.
    
    
    # By default, minimize
    if logical_not(exist('side','var')) or isempty(side):
        side='min'
# hessianextreme.m:65
    
    
    # If no initial guess was specified, prepare the empty one.
    if logical_not(exist('y0','var')):
        y0=matlabarray([])
# hessianextreme.m:70
    
    # Merge default solver options with potential user-specified options.
    # Set local defaults here
    localdefaults.verbosity = copy(0)
# hessianextreme.m:75
    localdefaults.solver = copy(trustregions)
# hessianextreme.m:76
    if logical_not(exist('options','var')) or isempty(options):
        options=struct()
# hessianextreme.m:78
    
    options=mergeOptions(localdefaults,options)
# hessianextreme.m:80
    
    if logical_not(exist('key','var')):
        if logical_not(exist('storedb','var')):
            storedb=StoreDB()
# hessianextreme.m:85
        key=storedb.getNewKey()
# hessianextreme.m:87
    
    
    # Convert the side into a sign.
    # Since Manopt minimizes, 'min' asks for no sign change.
    if 'min' == lower(side):
        sign=+ 1
# hessianextreme.m:94
    else:
        if 'max' == lower(side):
            sign=- 1
# hessianextreme.m:96
        else:
            error('The side should be either \'min\' or \'max\'.')
    
    # We define a manifold that is actually the unit sphere on the tangent
    # space to problem.M at x. A generalization would be to consider
    # Stiefel or Grassmann on the tangent space, but this would require
    # manipulating collections of tangent vectors, which in full generality
    # may be more complex (from a programming point of view).
    # Points are represented as tangent vectors of unit norm.
    # Tangent vectors are represented as tangent vectors orthogonal to the
    # root point, with respect to the Riemannian metric on the tangent
    # space.
    
    # M is the original manifold. x is a point on M.
    M=problem.M
# hessianextreme.m:112
    
    # tangent vector to M at x. This is a typical Riemannian submanifold of
    # a Euclidean space, hence it is easy to describe in terms of the tools
    # available for M.
    N=tangentspherefactory(M,x)
# hessianextreme.m:118
    
    # sure precomputable things are precomputed.
    if canGetGradient(problem):
        unused1,unused2=getCostGrad(problem,x,storedb,key,nargout=2)
# hessianextreme.m:123
    
    
    # This is the star operator of this party.
    hessian=lambda y=None: getHessian(problem,x,y,storedb,key)
# hessianextreme.m:127
    
    # problem on the sphere N.
    new_problem.M = copy(N)
# hessianextreme.m:131
    
    new_problem.cost = copy(cost)
# hessianextreme.m:135
    
@function
def cost(y=None,store=None,*args,**kwargs):
    varargin = cost.varargin
    nargin = cost.nargin

    store=prepare(y,store)
# hessianextreme.m:137
    f=dot(sign,store.f)
# hessianextreme.m:138
    return f,store
    
if __name__ == '__main__':
    pass
    
    new_problem.grad = copy(grad)
# hessianextreme.m:141
    
@function
def grad(y=None,store=None,*args,**kwargs):
    varargin = grad.varargin
    nargin = grad.nargin

    store=prepare(y,store)
# hessianextreme.m:143
    g=N.lincomb(y,dot(sign,2),store.Hy,dot(dot(sign,(- 2)),store.f),y)
# hessianextreme.m:144
    return g,store
    
if __name__ == '__main__':
    pass
    
    new_problem.hess = copy(hess)
# hessianextreme.m:147
    
@function
def hess(y=None,ydot=None,store=None,*args,**kwargs):
    varargin = hess.varargin
    nargin = hess.nargin

    store=prepare(y,store)
# hessianextreme.m:149
    Hydot=hessian[ydot]
# hessianextreme.m:150
    h=N.lincomb(y,dot(sign,2),Hydot,dot(dot(sign,(- 2)),store.f),ydot)
# hessianextreme.m:151
    h=N.proj(y,h)
# hessianextreme.m:152
    return h,store
    
if __name__ == '__main__':
    pass
    
    # This helper makes sure we do not duplicate Hessian computations.
    
@function
def prepare(y=None,store=None,*args,**kwargs):
    varargin = prepare.varargin
    nargin = prepare.nargin

    if logical_not(isfield(store,'ready')):
        Hy=hessian[y]
# hessianextreme.m:158
        store.f = copy(M.inner(x,y,Hy))
# hessianextreme.m:159
        store.Hy = copy(Hy)
# hessianextreme.m:160
        store.ready = copy(true)
# hessianextreme.m:161
    
    return store
    
if __name__ == '__main__':
    pass
    
    
    # Call a Manopt solver to solve the quadratic optimization problem on
    # the abstract sphere N.
    y,lambda_,info=manoptsolve(new_problem,y0,options,nargout=3)
# hessianextreme.m:167
    lambda_=dot(sign,lambda_)
# hessianextreme.m:168
    return store
    
if __name__ == '__main__':
    pass
    