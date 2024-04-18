import pandas as pd
import numpy as np

'''eand = function(coact,regDiscExp,multip=1)'''

def eand(coact,regDiscExp,multip=1):
    '''
    do.call(rbind,lapply(coact,function(co){
        if(co[1] == ""){
            return(rep(0,ncol(regDiscExp)))
        }else if(length(co) ==1){return(regDiscExp[co,])}
        n =length(co)
        x=apply(regDiscExp[co,],2,sum)
        y=x
        x[1:length(x)]=0
        x[which(y== - n )]=- multip
        x[which(y == n)] = multip
        return(x)
    }))
    '''
    def function_co(co):
        if co[1] == "":
            return np.zeros(regDiscExp.shape[1])
        elif len(co)==1:
            return regDiscExp[co, :]
        n=len(co)
        x = regDiscExp.iloc[co, :].sum(axis=0)
        y=x
        x[1:len(x)]=0
        x.loc[y[y==-n].index]=-multip
        x.loc[y[y==n].index]=multip
        return x

    concatened_df = pd.concat(map(lambda co: function_co(co),coact))
    return concatened_df

