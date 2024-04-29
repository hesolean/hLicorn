# librabies
import numpy as np

def eand(co_act, reg_disc_exp, multip=1) :
    def function_co(co) :
        if co['itemsets'] == "" :
            return np.zeros(len(reg_disc_exp.columns))
        elif len(co['itemsets']) == 1 :
            return reg_disc_exp.loc[co]
        else :
            n=len(co['itemsets'])
            y=reg_disc_exp.loc[co].apply(sum, axis=0)
            x=np.zeros(len(y))
            x[np.where(y == -n)[0]]=- multip
            x[np.where(y == n)[0]]=multip
            return x
    
    concatened_df=np.vstack([function_co(co) for co in co_act])
    return concatened_df