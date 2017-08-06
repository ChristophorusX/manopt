# Autogenerated with SMOP 
from smop.core import *
# getCost.m

    
@function
def getCost(problem=None,x=None,storedb=None,key=None,*args,**kwargs):
    varargin = getCost.varargin
    nargin = getCost.nargin

    # Computes the cost function at x.
    
    # function cost = getCost(problem, x)
# function cost = getCost(problem, x, storedb)
# function cost = getCost(problem, x, storedb, key)
    
    # Returns the value at x of the cost function described in the problem
# structure.
    
    # storedb is a StoreDB object, key is the StoreDB key to point x.
    
    # See also: canGetCost
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Dec. 30, 2012.
# Contributors: 
# Change log:
    
    #   April 3, 2015 (NB):
#       Works with the new StoreDB class system.
    
    # Allow omission of the key, and even of storedb.
    if logical_not(exist('key','var')):
        if logical_not(exist('storedb','var')):
            storedb=StoreDB()
# getCost.m:26
        key=storedb.getNewKey()
# getCost.m:28
    
    if isfield(problem,'cost'):
        ## Compute the cost function using cost.
        # Check whether this function wants to deal with storedb or not.
        if 1 == nargin(problem.cost):
            cost=problem.cost(x)
# getCost.m:38
        else:
            if 2 == nargin(problem.cost):
                # Obtain, pass along, and save the store for x.
                store=storedb.getWithShared(key)
# getCost.m:41
                cost,store=problem.cost(x,store,nargout=2)
# getCost.m:42
                storedb.setWithShared(store,key)
            else:
                if 3 == nargin(problem.cost):
                    # Pass along the whole storedb (by reference), with key.
                    cost=problem.cost(x,storedb,key)
# getCost.m:46
                else:
                    up=MException('manopt:getCost:badcost','cost should accept 1, 2 or 3 inputs.')
# getCost.m:48
                    throw(up)
    else:
        if isfield(problem,'costgrad'):
            ## Compute the cost function using costgrad.
            # Check whether this function wants to deal with storedb or not.
            if 1 == nargin(problem.costgrad):
                cost=problem.costgrad(x)
# getCost.m:59
            else:
                if 2 == nargin(problem.costgrad):
                    # Obtain, pass along, and save the store for x.
                    store=storedb.getWithShared(key)
# getCost.m:62
                    cost,grad,store=problem.costgrad(x,store,nargout=3)
# getCost.m:63
                    storedb.setWithShared(store,key)
                else:
                    if 3 == nargin(problem.costgrad):
                        # Pass along the whole storedb (by reference), with key.
                        cost=problem.costgrad(x,storedb,key)
# getCost.m:67
                    else:
                        up=MException('manopt:getCost:badcostgrad','costgrad should accept 1, 2 or 3 inputs.')
# getCost.m:69
                        throw(up)
        else:
            ## Abandon computing the cost function.
            up=MException('manopt:getCost:fail',cat('The problem description is not explicit enough to ','compute the cost.'))
# getCost.m:77
            throw(up)
    
    
    return cost
    
if __name__ == '__main__':
    pass
    