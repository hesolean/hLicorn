# librabies
import numpy as np

def eand(co_act, reg_disc_exp, multip=1) :
    concatened_df=[]
    for co in co_act :
        if co == "" :
            concatened_df.append(np.zeros(len(reg_disc_exp.columns)))
        elif len(co) == 1 :
            concatened_df.append(reg_disc_exp.loc[co])
        else :
            n=len(co)
            y=reg_disc_exp.loc[co].apply(sum, axis=0)
            x=np.zeros(len(y))
            x[np.where(y == -n)[0]]=- multip
            x[np.where(y == n)[0]]=multip
            concatened_df.append(x)
    concatened_df=np.vstack(concatened_df)
    return concatened_df