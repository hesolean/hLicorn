# librabies
import pandas as pd
import numpy as np
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori
import warnings
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

# secondary functions
from utils.discretizeExpressionData import discretizeExpressionData
from utils.oneGeneHLICORN import oneGeneHLICORN
    
def hLICORN(numericalExpression,Tflist,
            discreteExpression=None,GeneList=None,parallel="no",cluster=None,
            minGeneSupport=0.1,minCoregSupport = 0.1,maxCoreg=None,searchThresh=1/3,nGRN=100,verbose=False):
    # list of Tf
    tf_list = Tflist.iloc[:, 1].tolist()
    print(f"ft_list de taille ", len(tf_list))

    # determination of discrete values
    if discreteExpression==None:
        discreteExpression=discretizeExpressionData(numericalExpression)
    
    print(f"numericalExpression : {numericalExpression} de taille ", numericalExpression.shape)

    print(f"discreteExpression : {discreteExpression} de taille ", discreteExpression.shape)

    # list of row names of discreteExpression
    dis_row_names = None
    if not discreteExpression.empty:
        # access to names of the rows of discreteExpression to get a list
        dis_row_names = discreteExpression.index.tolist()
    print(f"dis_row_names de taille ", len(dis_row_names))

    # determination of GeneList
    num_row_names = None
    GeneList = None
    if GeneList==None:
        # access to names of the rows of numericalExpression to get a list
        num_row_names = numericalExpression.index.tolist()
        # get genes list witch are not in Tflist
        GeneList = list(set(num_row_names) - set(tf_list))
    print(f"GeneList de taille ", len(GeneList))

    # determination of maxCoreg
    if maxCoreg==None:
        maxCoreg=len(tf_list)
    print(f"maxCoreg : {maxCoreg}")

    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION BEFORE STARTING

    if  any(value not in [-1,0,1] for value in discreteExpression.stack().unique()):
        raise ValueError("Discrete expression data should only have values in {-1, 0, 1}")

    if len(num_row_names) > len(set(num_row_names)):
        raise ValueError("No gene duplicates are allowed in the row names.")

    # copy of the 2 lists to evaluate the condition without modify the original lists
    num_row_names_sort = num_row_names.copy()
    dis_row_names_sort = dis_row_names.copy()
    if numericalExpression.shape[0] != discreteExpression.shape[0] or num_row_names_sort.sort() != dis_row_names_sort.sort():
        raise ValueError("Discrete expression and continuous expression should have the same dimensions and the same rownames (gene/tf names)")

    if len(set(num_row_names).intersection(set(tf_list)))<=1:
        raise ValueError("At least 2 of the provided regulators/transcription factor (TFlist) should be in the rownames in the gene expression matrix")

    if len(set(num_row_names).intersection(set(GeneList)))==0:
        raise ValueError("The list of genes (GeneList) should be in the rownames in the gene expression matrix")

    parallel_options = ["multicore","no", "snow"]
    if parallel not in parallel_options:
        raise ValueError("The option of parallelism should be \"multicore\", \"no\" or \"snow\"")

    if verbose:
        print("Pre-process.")
    
    
    # select only genes and TF with ones or minus ones in at least minGeneSupport portion of samples
    genes_support = np.where(np.sum(np.abs(discreteExpression), axis=1) > (np.shape(numericalExpression)[1] * minGeneSupport))
    print(f"genes_support de taille ", genes_support[0].shape)

    discreteExpression=discreteExpression.iloc[genes_support]
    print(f"selected discreteExpression : {discreteExpression} de taille ", discreteExpression.shape)

    numericalExpression=numericalExpression.iloc[genes_support]
    print(f"selected numericalExpression de taille ", numericalExpression.shape)

    tf_list = list(set(num_row_names).intersection(tf_list))
    print(f"selected tf_list de taille ", len(tf_list))

    GeneList= list(set(num_row_names).intersection(GeneList))
    print(f"selected GeneList de taille ", len(GeneList))

    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION AFTER THE SELECTION OF THE GENES AND TF

    if len(tf_list)<5:
        raise ValueError("Less than 5 of the provided TF are suitable to infer a network. Either provide more TF, more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to select more but less variant TFs.")

    if len(GeneList)==0:
        raise ValueError("No genes were suitable to infer regulators. Either provide more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to allow the selection of more but less variant Genes.")
    
    
    # Get all the matrices and datasets needed (gene and tf expression, numerical or discrete)
    # If only one gene is given, R will automatically make a vector. The following make sure this does not happen.
    geneNumExp = numericalExpression.loc[GeneList]
    geneDiscExp= discreteExpression.loc[GeneList] 
    print("geneNumExp de taille ", geneNumExp.shape, " et geneDiscExp de taille ", geneDiscExp.shape)
    print(f"geneNumExp : {geneNumExp}")

    regNumExp= numericalExpression.loc[tf_list]
    regDiscExp= discreteExpression.loc[tf_list]
    print("regNumExp de taille ", regNumExp.shape, " et regDiscExp de taille ", regDiscExp.shape)
    print(f"regDiscExp : {regDiscExp}")


    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## TRANSFORMING ALL DISCRETE DATA INTO TRANSACTIONS
    # To run apriori, the discrete data must be binary. So, the discrete data is simply becoming two concatenated binary matrix
    # first n samples are positive expression values, then all negative values.

    posSamples = range(1,discreteExpression.shape[1]) # non utilisé

    negSamples= range(discreteExpression.shape[1]+1,discreteExpression.shape[1] *2) # non utilisé

    posRegDiscExp = regDiscExp.applymap(lambda x: 1 if x == -1 else x) # get a binary df with minus ones replaced by ones
    print(f"posRegDiscExp : {posRegDiscExp}")

    transRegBitData= posRegDiscExp.T
    print(f"transRegBitData : {transRegBitData} de type : ",type(transRegBitData))

    if verbose:
        print("Mining coregulator ...")
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## MINING FOR FREQUENT COREGULATORS
    # using apriori instead of eclat. testing may be required for possible speed improvement.
    '''
    miningFunction=apriori
    transitemfreq =suppressWarnings(miningFunction(transRegBitData,parameter=list(support = minGeneSupport/2,maxlen=1,target="frequent itemsets")
    ,control=list(verbose=FALSE))) -> on ne cherche que les singletons
    if(maxCoreg > 1){
        transitemfreq=c(transitemfreq,suppressWarnings(miningFunction(transRegBitData,parameter=list(support =minCoregSupport/2,minlen=2,maxlen=maxCoreg,target="closed frequent itemsets")
        ,control=list(verbose=FALSE))))
    } -> on ne cherche les "close" que dans un cas particulier
    coregs =as(slot(transitemfreq,"items"),"list")
    '''

    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        transitemfreq=apriori(transRegBitData, min_support=0.5, use_colnames=True, max_len=1, verbose=0, low_memory=True)
        # low_memory=True should only be used for large dataset if memory resources are limited 
    
        if maxCoreg > 1:
            result=apriori(transRegBitData, min_support=minCoregSupport/2, use_colnames=True, max_len=maxCoreg, verbose=0, low_memory=True)
            transitemfreq = pd.concat([transitemfreq, result])
        # To avoid duplicates in singletons
        coregs = set(transitemfreq['itemsets'].tolist())

    print(f"transitemfreq : {transitemfreq} de taille ", transitemfreq.shape)
    print(f"coregs : {coregs} de taille ", len(coregs))

    '''
    transitemfreq=c(transitemfreq,suppressWarnings(miningFunction(transRegBitData,parameter=list(support =minCoregSupport/2,minlen=2,maxlen=maxCoreg,target="closed frequent itemsets")
            ,control=list(verbose=FALSE))))
    J'ai un soucis avec le "closed frequent itemset" car normalement dans ce cas, il faut que maxlen soit à 0 ??
    minlen n'est pas un paramètre de apriori donc il faut le passer autrement ?

    '''


    '''
    if(verbose){
        message(paste("Learning a Co-Regulatory network for:\n",  length(GeneList)," target genes, ",length(TFlist)," regulators and a ",
        "total of coregulator sets ",length(coregs),"sets of potential co-regulators.\nSearch parameters :\n",
        "Maximum size of co-regulator sets : ",maxCoreg,"\nNumber of putative GRN per gene : ",
        nGRN,"\nMinimum number of differentially expressed samples to select a single gene : "
        ,minGeneSupport,"\nMinimum number of differentially expressed samples to select a set of co-regulator : "
        ,minCoregSupport,collapse=""))
        
        message("Mining GRN ...")
    }'''
    if verbose:
        print("Learning a Co-Regulatory network for:\n"+str(len(GeneList))
              +" target genes, "+str(len(tf_list))
              +" regulators and a total of coregulator sets "+str(len(coregs))
              +" sets of potential co-regulators.\nSearch parameters :\n"+"Maximum size of co-regulator sets : "+str(maxCoreg)
              +"\nNumber of putative GRN per gene : "+str(nGRN)+"\nMinimum number of differentially expressed samples to select a single gene : "+str(minGeneSupport)+
              "\nMinimum number of differentially expressed samples to select a set of co-regulator : "+str(minCoregSupport))
        
        print("Mining GRN ...")
    
    '''
    result=data.frame()
    gotNet=FALSE
    '''
    gotNet=False # on commence à un threathold élevé. Si on a pas de résultat, on réduit. C'est le gotNet qui dit si on continue ou non

    #just because it's easier toadd here 5% and remove it at the first line in the while loop, where it needs to be decrementale in case no GRNs are found
    '''searchThresh=  1/((1/searchThresh)-1)'''
    searchThresh=  1/((1/searchThresh)-1)
    print(f"searchThresh : {searchThresh}")
    # In very large datasets of very heterogeneous samples (such as the large collection of unrelated cell lines ...)
    # It is possible that no GRN can be fitted with stringent threshold (usually 50%) and that no GRN is found.
    # In case this happens, the threshold is decremented step by step and if no network is found at 10%, then none can be found ...
    '''
    while(searchThresh >= 0.05 & !gotNet )
    {
        searchThresh =  1/((1/searchThresh)+1)
        if(parallel =="multicore" & length(GeneList)>1 & getOption("mc.cores", 2L) > 1) -> joblib
        {
            result =mclapply(GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE -> a voir pour l'incorporer dans le oneGeneHLicorn
        }else if(parallel =="snow" & !is.null(cluster) & length(GeneList)>1){ -> spark
            result =parLapply(cluster,GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else if( length(GeneList)>1){ -> sequentiel
            result =lapply(GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else{  -> sequentiel
            result=oneGeneHLICORN(GeneList,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh ,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }
    }'''
    while searchThresh >= 0.05 and gotNet==False :
        print("Boucle while")
        # decrements the search threshold in case nothing is found
        #(can be the case for VERY large datasets for which it can be hard to find regulators with 50% of matching +1 and -1)
        searchThresh =  1/((1/searchThresh)+1)
        print(f"searchThresh : {searchThresh}")
        #running hlicorn for each gene in a multithread way if needed.
        
        # get available cares on the machine and define how much to use
        # available_cores_nb = multiprocessing.cpu_count()
        # using_processus = max(1, available_cores_nb - 1) 

        # on fait tourner joblib même si 1 coeur
        using_processus = 2
        print(f"using_processus : {using_processus}")

        def process_gene(gene, geneDiscExp, regDiscExp, coregs, transitemfreq, transRegBitData, searchThresh, geneNumExp, regNumExp, nGRN):
            print("process_gene")
            return oneGeneHLICORN(gene, geneDiscExp, regDiscExp, coregs, transitemfreq, transRegBitData, searchThresh, genexp=geneNumExp, regnexp=regNumExp, nresult=nGRN)

        if parallel =="multicore" and len(GeneList)>1 & using_processus > 1:
            print("multicore")
                        
            '''
            Version avec joblib
            '''            
            processes = [] # facultatif, uniquement pour récupérer la liste des processus qui se sont exécutés
            for _ in GeneList:
                p = multiprocessing.Process(target=process_gene,args=[GeneList])
                p.start()
                processes.append(p)
            
            for process in processes:
                process.join() # pour forcer le programme à attendre que tous les process soient finis avant de continuer
                
            # peu importe le nombre de coeurs disponibles, les taches seront distribuées en fonction des ressources disponibles.

            gotNet=True

        elif parallel =="snow" and cluster!=None & len(GeneList)>1:
            print("snow")
            with multiprocessing.Pool() as pool:
                results = [pool.apply_async(process_gene, (gene, geneDiscExp, regDiscExp, coregs, transitemfreq, transRegBitData, searchThresh, geneNumExp, regNumExp, nGRN)) for gene in GeneList]
                results = [res.get() for res in results]
                print(f"result : {results}")

            '''
            Pool() crée un pool de processus avec un nombre de processus par défaut (nombre de cœurs du CPU).
            pool.apply_async est utilisé pour exécuter de manière asynchrone la fonction process_gene pour chaque élément de GeneList.
            res.get() est utilisé pour obtenir les résultats de chaque tâche asynchrone.
            
            '''
            gotNet=True

        elif len(GeneList)>1:
            print("elif")
            # comprehension of the list to apply prosses_gene to each gene of the list
            results = [process_gene(gene, geneDiscExp, regDiscExp, coregs, transitemfreq, transRegBitData, searchThresh, geneNumExp, regNumExp, nGRN) for gene in GeneList]
            print(f"result : {results}")
            
            gotNet=True

        else:
            print("else")
            results=oneGeneHLICORN(GeneList,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            print(f"result : {results}")

            gotNet=True

    #if one of the above call to LICORN worked ... and if the result is a list (neither matrix nor dataframe)
    #  Merge the results into a data.Frame
    '''
    if(gotNet){
        #if needed, like when used in parallel, merge the results into a data.frame
        if(!is.data.frame(result) & !is.matrix(result)){result= data.frame(do.call(rbind,result))}
        #if LICORN actually did find some networks ... (meaning at least one GRN)
        if(ncol(result) >= 3 & nrow(result) >0){
            # Maybe LICORN did find somes nets, but not enough .. (for less then 5% of the genes)
            if(length(unique(result$Target)) < (0.05*length(GeneList))){
                gotNet=FALSE
            }
        }else{
            gotNet=FALSE
        }
        if(verbose){print(paste("got",nrow(result),"grn"))}
    }
    '''
    if gotNet:
        #if needed, like when used in parallel, merge the results into a data.frame
        if not isinstance(results, pd.DataFrame):
            '''& !is.matrix(result) : doit-on réellement tester si le résultat est une matrice puisqu'on est dans python ?'''
            results = pd.DataFrame(results)

        #if LICORN actually did find some networks ... (meaning at least one GRN)
        if results.shape[1] >= 3 and results.shape[0] >0:
            # Maybe LICORN did find somes nets, but not enough .. (for less then 5% of the genes)
            '''
            if(length(unique(result$Target)) < (0.05*length(GeneList))):
                gotNet=FALSE
            
            a-t-on une colonne Target dans le dataframe results ????'''
        else:
            gotNet=False
        
        if verbose:
            print("got"+str(results.shape[0])+"grn")
    
    '''
    if(verbose){
        message(paste("adjusted thresh:", searchThresh))
    }
    '''
    if verbose:
        print("adjusted thresh:"+str(searchThresh))
    
    
    #When done decrementing the threshold .... well if nothing was found maybe there is a probleme somewhere ...
    '''
    if(nrow(result) ==0 | ncol(result) <3){
        stop("Something went wrong. No GRN found.")
    }
    rownames(result)=NULL
    sigrns = coregnet(result)
    sigrns@inferenceParameters=list(minGeneSupport=minGeneSupport,maxCoreg=maxCoreg,minCoregSupport = minCoregSupport,searchThresh=searchThresh,nGRN=nGRN)
    return(sigrns)'''
    # if results.shape[0]==0 | results.shape[1] <3:
    #     raise ValueError("Something went wrong. No GRN found.")

    # results.iloc[:, 0] = None

    '''        FIN A ECLAIRCIR
        sigrns = coregnet(result)
    sigrns@inferenceParameters=list(minGeneSupport=minGeneSupport,maxCoreg=maxCoreg,minCoregSupport = minCoregSupport,searchThresh=searchThresh,nGRN=nGRN)
        return(sigrns)'''
        