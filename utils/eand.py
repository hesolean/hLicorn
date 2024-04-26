# librabies
import numpy as np

def eand(coact, regDiscExp, multip=1):
    def function_co(co):
        print("type de coact : ",type(coact))
        if co['itemsets'] == "" :
            return np.zeros(len(regDiscExp.columns))
        elif len(co['itemsets']) == 1:
            return regDiscExp.loc[co]
        else:
            n = len(co['itemsets'])
            y = regDiscExp.loc[co].apply(sum, axis=0)
            x = np.zeros(len(y))
            x[np.where(y == -n)[0]] = - multip
            x[np.where(y == n)[0]] = multip
            return x
    
    concatened_df = np.vstack([function_co(co) for co in coact])
    return concatened_df