# librabies
import numpy as np
from numpy import intersect1d
import pandas as pd
from sklearn.linear_model import LinearRegression

# secondary functions
from utils.eand import eand
from utils.C_call import C_call

# function that return new frequent itemsets according to the studied gene
def support(items, test_df) :
    # print to check
    print("support")

    freq_items=pd.DataFrame()
    supports=[]

    # map on the itemset already find in hLicorn
    for itemset in items.itemsets :
        freq=test_df.shape[0]
        for item in itemset :
            sup=0
            if item in test_df.columns :	sup = test_df[test_df[item] == True][item].count()
            if sup<freq :	freq=sup
        # print(str(itemset), " : " , str(freq))
        supports.append(freq / test_df.shape[0])
        # print("supports ", supports)

    # create new frequent itemset
    freq_items["supports"]=supports
    freq_items["items"]=items.itemsets

    return freq_items

def oneGeneHLICORN(g, gene_disc_exp, reg_disc_exp, co_regs, trans_item_freq, trans_reg_bit_data, search_thresh, reg_num_exp, gene_num_exp, nresult) :
    '''try to run as R code ie with index'''
    # shift=gene_disc_exp.shape[1]

    # sample index with the target gene at 1 or -1
    # pos=np.where(gene_disc_exp.loc[g, :] == 1)[0]
    # neg=np.where(gene_disc_exp.loc[g, :] == -1)[0] + shift

    # co_act_fi=support(trans_item_freq, trans_reg_bit_data.iloc[:, pos + neg])
    # co_rep_fi=support(trans_item_freq, trans_reg_bit_data.iloc[:, neg + pos])

    '''Alternative with column names instead of index'''
    # retrieve the column names corresponding to the values of 1 in gene_disc_exp
    pos_columns=gene_disc_exp.columns[gene_disc_exp.loc[g] == 1]
    neg_columns=gene_disc_exp.columns[gene_disc_exp.loc[g] == 1]

    # add "pos_" befor the name of the columns to have no conflict of index
    pos_columns_with_prefix = [f"pos_{col}" for col in pos_columns]
    neg_columns_with_prefix=[f"neg_{col}" for col in neg_columns]

    # print to check
    print("pos", pos_columns_with_prefix)
    print("neg ", neg_columns_with_prefix)
    print("trans_reg_bit_data.iloc[:, pos + neg]", trans_reg_bit_data.loc[pos_columns_with_prefix + neg_columns_with_prefix, :])
    # select all the coregulators with a support of 50% minimum only in the samples with the target gene at ones or minus ones
    # indices for which the threshold is reached, then we select the elements
    
    co_act_fi=support(trans_item_freq, trans_reg_bit_data.loc[pos_columns_with_prefix + neg_columns_with_prefix, :])
    co_act = [item for item, support in zip(co_regs, co_act_fi.supports) if support >= search_thresh]

    co_rep_fi=support(trans_item_freq, trans_reg_bit_data.loc[neg_columns_with_prefix + pos_columns_with_prefix, :])
    co_rep=[item for item, support in zip(co_regs, co_rep_fi.supports) if support >= search_thresh]

    # to have unique coregulators and a single vector of coreg (not a list)
    # extraction of the names in each frozenset
    co_act_names=sorted(list(set([item for sublist in co_act for item in sublist])))
    co_rep_names=sorted(list(set([item for sublist in co_rep for item in sublist])))
    
    # print to check
    print("co_act_names : ", co_act_names)

    # add empty coregulators to have the possibility to only have ativators or inhibitors
    co_rep.append("")
    co_act.append("")

    # merge expression of coregulator and corepressor
    co_act_exp=eand(co_act, reg_disc_exp)
    co_rep_exp=eand(co_rep, reg_disc_exp)
    
    # active inhibitor has a stronger impact than naything else in licorn:
    co_rep_exp=co_rep_exp * 2

    # print to check
    print("co_rep_exp : ", co_rep_exp)

    C_call(co_act_exp, co_act, co_rep_exp, co_rep, g, gene_disc_exp)

    # bad index will store all bad comparisons (could be done before computing .. right?)
    # then, no intersection between act and rep (includes both empty

    '''TRANSLATION OF THE R CODE : PART TO VERIFY ACCORDING TO THE NEEDS'''
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