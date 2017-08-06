# Autogenerated with SMOP 
from smop.core import *
# applyStatsfun.m

    
@function
def applyStatsfun(problem=None,x=None,storedb=None,key=None,options=None,stats=None,*args,**kwargs):
    varargin = applyStatsfun.varargin
    nargin = applyStatsfun.nargin

    # Apply the statsfun function to a stats structure (for solvers).
    
    # function stats = applyStatsfun(problem, x, storedb, key, options, stats)
    
    # Applies the options.statsfun user supplied function (if it was provided)
# to the stats structure, and returns the (possibly) modified stats
# structure.
    
    # storedb is a StoreDB object, key is the StoreDB key to point x.
    
    # Note: if statsfun accepts a store structure as input, this structure can
# be read but not modified (modifications will be lost) ; the store
# structure will contain the store.shared field.
    
    # See also:
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, April 3, 2013.
# Contributors: 
# Change log:
    
    #   April 3, 2015 (NB):
#       Works with the new StoreDB class system.
    
    if isfield(options,'statsfun'):
        if 3 == nargin(options.statsfun):
            stats=options.statsfun(problem,x,stats)
# applyStatsfun.m:30
        else:
            if 4 == nargin(options.statsfun):
                # Obtain, pass along, and save the store for x.
                # get/setWithShared must come in pairs.
                store=storedb.getWithShared(key)
# applyStatsfun.m:34
                stats=options.statsfun(problem,x,stats,store)
# applyStatsfun.m:35
                storedb.setWithShared(store,key)
            else:
                warning('manopt:statsfun','statsfun unused: wrong number of inputs')
    
    return stats
    
if __name__ == '__main__':
    pass
    