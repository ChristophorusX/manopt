# Autogenerated with SMOP 
from smop.core import *
# tangentspacefactory.m

    
@function
def tangentspacefactory(M=None,x=None,*args,**kwargs):
    varargin = tangentspacefactory.varargin
    nargin = tangentspacefactory.nargin

    # Returns a manifold structure representing the tangent space to M at x.
    
    # N = tangentspacefactory(M, x)
    
    # N defines a (linear) manifold that is the tangent space to M at x. Points
# are represented as tangent vectors to M at x. Tangent vectors are also
# represented as tangent vectors to M at x.
    
    # This is chiefly useful to solve optimization problems involving tangent
# vectors to M at x, which notably comes up when solving linear systems
# involving, for example, the Hessian of the cost on M at x. The Riemannian
# (actually, Euclidean) structure on N is that of the tangent space to M,
# that is, the inner product is inherited.
    
    # See also: preconhessiansolve
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, April 9, 2015.
# Contributors: 
# Change log:
    
    # N is the manifold we build. y will be a point on N, thus also a
    # tangent vector to M at x. This is a typical Euclidean space, hence it
    # will be easy to describe in terms of the tools available for M.
    N=struct()
# tangentspacefactory.m:26
    
    # N at y is the tangent space to M at x, thus u, u1 and u2 are also
    # tangent vectors to M at x.
    
    N.dim = copy(lambda : M.dim())
# tangentspacefactory.m:32
    N.inner = copy(lambda y=None,u1=None,u2=None: M.inner(x,u1,u2))
# tangentspacefactory.m:33
    N.norm = copy(lambda y=None,u=None: M.norm(x,u))
# tangentspacefactory.m:34
    N.proj = copy(M.proj)
# tangentspacefactory.m:35
    N.typicaldist = copy(lambda : N.dim())
# tangentspacefactory.m:36
    N.tangent = copy(lambda y=None,u=None: u)
# tangentspacefactory.m:37
    N.egrad2rgrad = copy(lambda x=None,g=None: g)
# tangentspacefactory.m:38
    N.ehess2rhess = copy(lambda x=None,eg=None,eh=None,d=None: eh)
# tangentspacefactory.m:39
    N.exp = copy(exponential)
# tangentspacefactory.m:40
    N.retr = copy(exponential)
# tangentspacefactory.m:41
    N.log = copy(lambda y1=None,y2=None: M.lincomb(x,1,y2,- 1,y1))
# tangentspacefactory.m:42
    N.pairmean = copy(lambda y1=None,y2=None: M.lincomb(x,0.5,y1,0.5,y2))
# tangentspacefactory.m:43
    N.rand = copy(lambda : M.randvec(x))
# tangentspacefactory.m:44
    N.randvec = copy(lambda y=None: M.randvec(x))
# tangentspacefactory.m:45
    N.zerovec = copy(M.zerovec)
# tangentspacefactory.m:46
    N.lincomb = copy(M.lincomb)
# tangentspacefactory.m:47
    N.transp = copy(lambda y1=None,y2=None,u=None: u)
# tangentspacefactory.m:48
    N.hash = copy(lambda y=None: cat('z',hashmd5(M.vec(x,y))))
# tangentspacefactory.m:49
    
    
@function
def exponential(y=None,u=None,t=None,*args,**kwargs):
    varargin = exponential.varargin
    nargin = exponential.nargin

    if nargin == 2:
        t=1
# tangentspacefactory.m:54
    
    yy=M.lincomb(x,1,y,t,u)
# tangentspacefactory.m:56
    return yy
    
if __name__ == '__main__':
    pass
    
    
    return yy
    
if __name__ == '__main__':
    pass
    