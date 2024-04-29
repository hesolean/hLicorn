# librabies
import numpy as np
from numpy import intersect1d
import pandas as pd
from sklearn.linear_model import LinearRegression

# secondary functions
from utils.eand import eand
from utils.C_call import C_call


def oneGeneHLICORN(g, gene_disc_exp, reg_disc_exp, co_regs, trans_reg_bit_data, search_thresh, gen_exp, nresult) :
    print("oneGene")
    shift=gene_disc_exp.shape[1]

    #sample index with the target gene at 1 or -1
    pos=np.where(gene_disc_exp.loc[g, :] == 1)[0]
    neg=np.where(gene_disc_exp.loc[g, :] == -1)[0] + shift

    # function that return new frequent itemset according to the studied gene
    def support(items, test_df) :
        freq_items=pd.DataFrame()
        supports=[]
        for itemset in items.itemsets :
            freq=test_df.shape[0]
            for item in itemset :
                sup=0
                if item in test_df.columns :	sup = test_df[test_df[item] == True][item].count()
                if sup<freq :	freq=sup
            print(str(itemset), " : " , str(freq))
            supports.append(freq / test_df.shape[0])

        freq_items["support"]=supports
        freq_items["items"]=items.itemsets
	
        return freq_items

    # select all the coregulators with a support of 50% minimum only in the samples with the target gene at ones or minus ones
    # indices for which the threshold is reached, then we select the elements
    co_act_fi=support(trans_reg_bit_data.iloc[:, pos + neg], g)
    co_act=[item for item in co_regs if co_act_fi >= search_thresh]
    co_rep_fi=support(trans_reg_bit_data.iloc[:, neg + pos], g)
    co_rep=[item for item in co_regs if co_rep_fi >= search_thresh]

    # add empty coregulators to have the possibility to only have ativators or inhibitors
    co_rep.loc[len(co_rep)]=""
    co_act.loc[len(co_act)]=""

    # to have unique coregulators and a single vector of coreg (not a list)
    co_act_names=sorted(set(co_act['itemsets'].split()))
    co_rep_names=sorted(set(co_rep['itemsets'].split()))
    co_act=co_act.drop_duplicates()
    co_rep=co_rep.drop_duplicates()

    # merge expression of coregulator and corepressor
    # print(f"reg_disc_exp : {reg_disc_exp}")
    co_act_exp=eand(co_act, reg_disc_exp)
    co_rep_exp=eand(co_rep, reg_disc_exp)
    
    # active inhibitor has a stronger impact than naything else in licorn:
    co_rep_exp.replace(1, 2, inplace=True)

    C_call(co_act_exp, co_act, co_rep_exp, co_rep, g, gene_disc_exp)

    # bad index will store all bad comparisons (could be done before computing .. right?)
    # then, no intersection between act and rep (includes both empty

    # empilement de x[8] et x[9] mais je ne sais pas à quoi correspondent les x[8], x[9]
    combined_matrix=np.column_stack((x[[8]], x[[9]]))

    good_index=[]

    # Boucle for pour appliquer intersect() à chaque ligne
    for i, column in enumerate(combined_matrix) :
        vector1=co_act[[column[0]]]
        vector2=co_rep[[column[1]]]
        intersection_length=len(intersect1d(vector1, vector2))
        if intersection_length == 0 :
            good_index.append(i)

    sel_act=[co_act_names[x[[8]]] for i in good_index]
    sel_rep=[co_rep_names[x[[9]]] for i in good_index]

    # all empty set of coregulators are set to NA
    sel_act[pd.isempty(sel_act)]=pd.NA
    sel_rep[pd.isempty(sel_rep)]=pd.NA

    mae=[x[[7]] for i in good_index]
    if not pd.isna(nresult) :
        # get 100 first ranks, if ties, might get more ...
        best_index=mae[mae.rank(method='min') <= nresult].index
        
    else :
        best_index=mae[mae.rank(method='min') == 1].index

    # new df with targets and co_regs
    GRN = pd.DataFrame({
        "Target": [g[i] for i in best_index for _ in range(len(g))],
        "co_act": [sel_act[i] for i in best_index],
        "co_rep": [sel_rep[i] for i in best_index]
    })

    # if no grn are found return NULL
    if GRN.shape == 0 :
        return None

    if isinstance(GRN, np.ndarray) or isinstance(GRN, pd.DataFrame) :
        linear_models=np.apply_along_axis(
            lambda row : LinearRegression().fit(row.reshape(-1, 1), gen_exp.reshape(-1, 1)).coef_[0], 
            axis=1, 
            arr=GRN
        )
    else :
        # if GRN is not matrix nor dataframe
        linear_models=LinearRegression().fit(GRN.reshape(-1, 1), gen_exp.reshape(-1, 1)).coef_[0]

    num_scores=pd.DataFrame({"numscores": linear_models})
    num_scores.iloc[:, 2:6]=num_scores.iloc[:, 2:6].apply(pd.to_numeric)
    num_scores.columns=["Target", "co_act", "co_rep"]

    result=pd.concat([GRN, num_scores], axis=1, ignore_index=True)

    # type conversion in string
    result=result.astype(str)

    # rename the columns
    result.columns=["Target", "co_act", "co_rep", "num_scores1", "num_scores2", "num_scores3", "num_scores4"]
    print("fin oneGene")
    return  result