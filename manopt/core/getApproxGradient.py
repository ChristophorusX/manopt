# Autogenerated with SMOP 
from smop.core import *
# getApproxGradient.m

    
@function
def getApproxGradient(problem=None,x=None,storedb=None,key=None,*args,**kwargs):
    varargin = getApproxGradient.varargin
    nargin = getApproxGradient.nargin

    # Computes an approximation of the gradient of the cost function at x.
    
    # function approxgrad = getApproxGradient(problem, x)
# function approxgrad = getApproxGradient(problem, x, storedb)
# function approxgrad = getApproxGradient(problem, x, storedb, key)
    
    # Returns an approximation of the gradient at x for the cost function
# described in the problem structure.
    
    # storedb is a StoreDB object, key is the StoreDB key to point x.
    
    # If no approximate gradient was provided, this call is redirected to
# getGradientFD.
# 
# See also: getGradientFD canGetApproxGradient
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Nov. 1, 2016.
# Contributors: 
# Change log:
    
    # Allow omission of the key, and even of storedb.
    if logical_not(exist('key','var')):
        if logical_not(exist('storedb','var')):
            storedb=StoreDB()
# getApproxGradient.m:26
        key=storedb.getNewKey()
# getApproxGradient.m:28
    
    if isfield(problem,'approxgrad'):
        ## Compute the approximate gradient using approxgrad.
        # Check whether this function wants to deal with storedb or not.
        if 1 == nargin(problem.approxgrad):
            approxgrad=problem.approxgrad(x)
# getApproxGradient.m:38
        else:
            if 2 == nargin(problem.approxgrad):
                # Obtain, pass along, and save the store for x.
                store=storedb.getWithShared(key)
# getApproxGradient.m:41
                approxgrad,store=problem.approxgrad(x,store,nargout=2)
# getApproxGradient.m:42
                storedb.setWithShared(store,key)
            else:
                if 3 == nargin(problem.approxgrad):
                    # Pass along the whole storedb (by reference), with key.
                    approxgrad=problem.approxgrad(x,storedb,key)
# getApproxGradient.m:46
                else:
                    up=MException('manopt:getApproxGradient:badapproxgrad','approxgrad should accept 1, 2 or 3 inputs.')
# getApproxGradient.m:48
                    throw(up)
    else:
        ## Try to fall back to a standard FD approximation.
        approxgrad=getGradientFD(problem,x,storedb,key)
# getApproxGradient.m:56
    
    
    return approxgrad
    
if __name__ == '__main__':
    pass
    