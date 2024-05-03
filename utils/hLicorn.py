# librabies
import pandas as pd
import numpy as np
from mlxtend.frequent_patterns import apriori
import warnings
import multiprocessing

# secondary functions
from utils.discretize_expression_data import discretize_expression_data
from utils.oneGeneHLICORN import oneGeneHLICORN

def hLICORN(numerical_expression, tf_list,
            discrete_expression=None, gene_list=None, parallel="no", cluster=None,
            min_gene_support=0.1, min_coreg_support = 0.1, max_coreg=None, search_thresh=1/3, nGRN=100, verbose=False) :
    # list of Tf
    tf_list=tf_list.iloc[:, 1].tolist()

    # determination of discrete values
    if discrete_expression == None :
        discrete_expression=discretize_expression_data(numerical_expression)

    # list of row names of discrete_expression
    dis_row_names=None
    if not discrete_expression.empty :
        # access to names of the rows of discrete_expression to get a list
        dis_row_names=discrete_expression.index.tolist()

    # determination of gene_list
    num_row_names=None
    if gene_list == None :
        # access to names of the rows of numerical_expression to get a list
        num_row_names=numerical_expression.index.tolist()
        # get genes list witch are not in tf_list
        gene_list=list(set(num_row_names) - set(tf_list))

    # determination of max_coreg if not define
    if max_coreg == None :
        max_coreg=len(tf_list)

    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION BEFORE STARTING

    if any(value not in [-1,0,1] for value in discrete_expression.stack().unique()):
        raise ValueError("Discrete expression data should only have values in {-1, 0, 1}")

    if len(num_row_names) > len(set(num_row_names)) :
        raise ValueError("No gene duplicates are allowed in the row names.")

    # copy of the 2 lists to evaluate the condition without modify the original lists
    num_row_names_sort=num_row_names.copy()
    dis_row_names_sort=dis_row_names.copy()
    if numerical_expression.shape[0] != discrete_expression.shape[0] or num_row_names_sort.sort() != dis_row_names_sort.sort() :
        raise ValueError("Discrete expression and continuous expression should have the same dimensions and the same rownames (gene/tf names)")

    if len(set(num_row_names).intersection(set(tf_list))) <= 1 :
        raise ValueError("At least 2 of the provided regulators/transcription factor (tf_list) should be in the rownames in the gene expression matrix")

    if len(set(num_row_names).intersection(set(gene_list))) == 0 :
        raise ValueError("The list of genes (gene_list) should be in the rownames in the gene expression matrix")

    # define options for the parallelization of the program
    parallel_options=["multicore", "no", "spark"]
    if parallel not in parallel_options :
        raise ValueError("The option of parallelism should be \"multicore\", \"no\" or \"spark\"")

    if verbose :
        print("Pre-process.")
    
    # select only genes and TF with ones or minus ones in at least min_gene_support portion of samples
    genes_support=np.where(np.sum(np.abs(discrete_expression), axis=1) > (np.shape(numerical_expression)[1] * min_gene_support))
    discrete_expression=discrete_expression.iloc[genes_support]
    numerical_expression=numerical_expression.iloc[genes_support]
    tf_list=list(set(num_row_names).intersection(tf_list))
    gene_list=list(set(num_row_names).intersection(gene_list))

    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION AFTER THE SELECTION OF THE GENES AND TF

    if len(tf_list) < 5:
        raise ValueError("Less than 5 of the provided TF are suitable to infer a network. Either provide more TF, more variations in the discrete dataset (more 1 or -1) or decrease the min_gene_support parameter to select more but less variant TFs.")

    if len(gene_list) == 0 :
        raise ValueError("No genes were suitable to infer regulators. Either provide more variations in the discrete dataset (more 1 or -1) or decrease the min_gene_support parameter to allow the selection of more but less variant Genes.")
    
    
    # Get all the matrix and datasets needed (gene and tf expression, numerical or discrete)
    # If only one gene is given, R will automatically make a vector. The following make sure this does not happen.
    gene_num_exp=numerical_expression.loc[gene_list]
    gene_disc_exp=discrete_expression.loc[gene_list] 

    reg_num_exp=numerical_expression.loc[tf_list]
    reg_disc_exp=discrete_expression.loc[tf_list]

    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## TRANSFORMING ALL DISCRETE DATA INTO TRANSACTIONS
    # To run apriori, the discrete data must be binary. So, the discrete data is simply becoming two concatenated binary matrix
    # first n samples are positive expression values, then all negative values.

    pos_samples=range(1,discrete_expression.shape[1]) # not used
    neg_samples=range(discrete_expression.shape[1]+1,discrete_expression.shape[1] *2) # not used

    # tranform minus ones in zeros and rename columns to get positives set
    pos_reg_disc_exp=reg_disc_exp.map(lambda x: 0 if x == -1 else x)
    pos_reg_disc_exp.columns = [f"pos_{col}" for col in pos_reg_disc_exp.columns]

    # transform ones in zeros and minus ones in ones for the apriori and rename columns to get negatives set
    neg_reg_disc_exp=reg_disc_exp.map(lambda x: 0 if x == 1 else 1 if x == -1 else x)
    neg_reg_disc_exp.columns = [f"neg_{col}" for col in neg_reg_disc_exp.columns]

    # concat the 2 dataframes in one with positive values then negative values
    global_reg_disc_exp=pd.concat([pos_reg_disc_exp, neg_reg_disc_exp], axis=1)

    # transpose for apriori function
    trans_reg_bit_data=global_reg_disc_exp.T

    # print to check
    print(f"trans_reg_bit_data :  {trans_reg_bit_data}")
    
    if verbose :
        print("Mining coregulator ...")
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## MINING FOR FREQUENT COREGULATORS
    # using apriori instead of eclat. testing may be required for possible speed improvement.
    with warnings.catch_warnings():
        # print to check
        print("Calcul des itemset frequents")
        warnings.simplefilter("ignore") 
        trans_item_freq=apriori(trans_reg_bit_data, min_support=0.9, use_colnames=True, max_len=1, verbose=0)
    
        if max_coreg > 1 :
            result=apriori(trans_reg_bit_data, min_support=min_coreg_support/2, use_colnames=True, max_len=max_coreg, verbose=0)
            trans_item_freq=pd.concat([trans_item_freq, result])
        
        # to avoid duplicates in singletons
        co_regs=set(trans_item_freq['itemsets'].tolist())

    if verbose :
        print("Learning a Co-Regulatory network for:\n" + str(len(gene_list))
              + " target genes, " + str(len(tf_list))
              + " regulators and a total of coregulator sets " + str(len(co_regs))
              + " sets of potential co-regulators.\nSearch parameters :\n" + "Maximum size of co-regulator sets : " + str(max_coreg)
              + "\nNumber of putative GRN per gene : " + str(nGRN) + "\nMinimum number of differentially expressed samples to select a single gene : " + str(min_gene_support)
              + "\nMinimum number of differentially expressed samples to select a set of co-regulator : " + str(min_coreg_support))
        
        print("Mining GRN ...")

    # initialization of boolean value
    got_net=False

    #just because it's easier to add here 5% and remove it at the first line in the while loop, where it needs to be decrementale in case no GRNs are found
    search_thresh=1 / ((1 / search_thresh) - 1)

    # In very large datasets of very heterogeneous samples (such as the large collection of unrelated cell lines ...)
    # It is possible that no GRN can be fitted with stringent threshold (usually 50%) and that no GRN is found.
    # In case this happens, the threshold is decremented step by step and if no network is found at 10%, then none can be found ...
    while search_thresh >= 0.05 and got_net == False :

        # decrements the search threshold in case nothing is found
        #(can be the case for VERY large datasets for which it can be hard to find regulators with 50% of matching +1 and -1)
        search_thresh=1 / ((1 / search_thresh) + 1)

        #running hlicorn for each gene in a multithread way if needed.
        
        # get available cares on the machine and define how much to use
        # with joblib, no matter number of processus, the distribution is automatic
        # available_cores_nb = multiprocessing.cpu_count()
        # using_processus = max(1, available_cores_nb - 1) 
        
        # fix using_processus to run even with less than 2 cores
        using_processus=2

        if parallel == "multicore" and len(gene_list) > 1 and using_processus > 1 :
            # print to check
            print("multicore")
            # joblib to code
            got_net=True

        elif parallel == "spark" :
            # print to check
            print("spark")
            # spark to code
            got_net=True

        elif len(gene_list) > 1 :
            # print to check
            print("3eme")
            results=[]

            # comprehension of the list to apply oneGeneHLICORN to each gene of the list
            for gene in gene_list :
                # print to check
                print(f"gene : {gene}")
                results.append(oneGeneHLICORN(gene, gene_disc_exp, reg_disc_exp, co_regs, trans_item_freq, trans_reg_bit_data, search_thresh, reg_num_exp, gene_num_exp, nresult=nGRN))
            
            got_net=True

        else :
            # print to check
            print("dernier")

            # in this case, gene_list only have one element
            gene=gene_list[0]
            results=oneGeneHLICORN(gene, gene_disc_exp, reg_disc_exp, co_regs, trans_item_freq, trans_reg_bit_data, search_thresh, reg_num_exp, gene_num_exp, nresult=nGRN)

            got_net=True


    '''        FIN A VERIFIER      '''


    # if one of the above call to LICORN worked ... and if the result is a list (neither matrix nor dataframe)
    # merge the results into a dataframe
    if got_net:
        # print to check
        print("got_net")

        # if needed, like when used in parallel, merge the results into a dataframe
        if not isinstance(results, pd.DataFrame) :
            results=pd.DataFrame(results)

        else :
            got_net=False

        # if LICORN actually did find some networks ... (meaning at least one GRN)
        if results.shape[1] >= 3 and results.shape[0] > 0 :
            # Maybe LICORN did find somes nets, but not enough .. (for less then 5% of the genes)
            if len(result['Target'].unique().tolist()) < (0.05 * len(gene_list)) :
                got_net=False
        else :
            got_net=False

        if verbose :
            print("got" + str(results.shape[0]) + "grn")
    
    if verbose:
        print("adjusted thresh:" + str(search_thresh))
    
    #When done decrementing the threshold .... well if nothing was found maybe there is a probleme somewhere ...
    if results.shape[0]==0 | results.shape[1] <3:
        raise ValueError("Something went wrong. No GRN found.")

    results.iloc[:, 0]=None
    results.index=None
    # sigrns=coregnet(results)
    # sigrns@inferenceParameters=list(minGeneSupport=minGeneSupport,maxCoreg=maxCoreg,minCoregSupport = minCoregSupport,searchThresh=searchThresh,nGRN=nGRN)
    # return sigrns
        