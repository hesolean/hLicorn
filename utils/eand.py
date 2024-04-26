# librabies
import pandas as pd
import numpy as np

def eand(coact, regDiscExp, multip=1):
    def function_co(co):
        print("type de coact : ",type(coact))
        if co[0] == False :
            return np.zeros(len(regDiscExp.columns))
        elif len(co) == 1:
            return regDiscExp.loc[co]
        else:
            y = regDiscExp.loc[co].apply(sum, axis=0)
            if y == -len(co):
                return -1 * multip
            elif y == len(co):
                return 1 * multip
            else:
                return y
    
    concatened_df = np.vstack([function_co(co) for co in coact])
    return concatened_df

