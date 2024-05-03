# librabies
import numpy as np

def eand(co_act, reg_disc_exp, multip=1) :
    # create df
    concatened_df=[]

    # map on co_act / co_rep df
    for co in co_act :
        if co == "" :
            # generate a linear array with zeros
            concatened_df.append(np.zeros(len(reg_disc_exp.columns)))

        elif len(co) == 1 :
            # get the co_act / co_rep as it is
            concatened_df.append(reg_disc_exp.loc[co])

        else :
            # generate a linear array with ones and minus ones based on the sum of identical co_act and co_rep
            n=len(co)
            y=reg_disc_exp.loc[co].apply(sum, axis=0)
            x=np.zeros(len(y))
            x[np.where(y == -n)[0]]= - multip
            x[np.where(y == n)[0]]=multip
            concatened_df.append(x)

    # vertical stacking of the df results
    concatened_df=np.vstack(concatened_df)
    
    return concatened_df