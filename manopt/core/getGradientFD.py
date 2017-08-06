# Autogenerated with SMOP 
from smop.core import *
# getGradientFD.m

    
@function
def getGradientFD(problem=None,x=None,storedb=None,key=None,*args,**kwargs):
    varargin = getGradientFD.varargin
    nargin = getGradientFD.nargin

    # Computes an approx. of the gradient w/ finite differences of the cost.
    
    # function gradfd = getGradientFD(problem, x)
# function gradfd = getGradientFD(problem, x, storedb)
# function gradfd = getGradientFD(problem, x, storedb, key)
    
    # Returns a finite difference approximation of the gradient at x for
# the cost function described in the problem structure. The finite
# difference is based on M.dim()+1 computations of the cost.
    
    # storedb is a StoreDB object, key is the StoreDB key to point x.
    
    # If the cost cannot be computed, an exception is thrown.
    
    # See also: approxgradientFD
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Nov. 1, 2016.
# Contributors: 
# Change log:
    
    # Allow omission of the key, and even of storedb.
    if logical_not(exist('key','var')):
        if logical_not(exist('storedb','var')):
            storedb=StoreDB()
# getGradientFD.m:26
        key=storedb.getNewKey()
# getGradientFD.m:28
    
    # This gradient approximation is based on the cost:
    # check availability.
    if logical_not(canGetCost(problem)):
        up=MException('manopt:getGradientFD:nocost','getGradientFD requires the cost to be computable.')
# getGradientFD.m:34
        throw(up)
    
    
    
    # Default parameters. See approxgradientFD for explicit user access to
    # these parameters.
    stepsize=2 ** - 23
# getGradientFD.m:42
    subspacedim=matlabarray([])
# getGradientFD.m:43
    
    fx=getCost(problem,x,storedb,key)
# getGradientFD.m:47
    
    # thereof. The default is a full subspace. If a strict subspace is
    # picked, the returned vector approximates the orthogonal projection of
    # the gradient to that subspace.
    B=tangentorthobasis(problem.M,x,subspacedim)
# getGradientFD.m:53
    
    # along each direction in the basis B.
    df=zeros(size(B))
# getGradientFD.m:57
    for k in arange(1,numel(B)).reshape(-1):
        # Move in the B{k} direction
        xk=problem.M.retr(x,B[k],stepsize)
# getGradientFD.m:60
        fxk=getCost(problem,xk,storedb)
# getGradientFD.m:62
        df[k]=(fxk - fx) / stepsize
# getGradientFD.m:64
    
    
    # Build the gradient approximation.
    gradfd=lincomb(problem.M,x,B,df)
# getGradientFD.m:68
    return gradfd
    
if __name__ == '__main__':
    pass
    