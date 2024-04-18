import pandas as pd
import numpy as np
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori
import warnings
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from utils.discretizeExpressionData import discretizeExpressionData
from utils.oneGeneHLICORN import oneGeneHLICORN

'''
hLICORN=function( numericalExpression,discreteExpression=discretizeExpressionData(numericalExpression)
, TFlist, GeneList=setdiff(rownames(numericalExpression),TFlist),parallel = c("multicore","no", "snow"),cluster=NULL,
minGeneSupport=0.1,minCoregSupport = 0.1,maxCoreg=length(TFlist),
searchThresh=1/3,nGRN=100,verbose=FALSE)
'''    
def hLICORN(numericalExpression,Tflist,
            discreteExpression=None,GeneList=None,parallel="no",cluster=None,
            minGeneSupport=0.1,minCoregSupport = 0.1,maxCoreg=None,searchThresh=1/3,nGRN=100,verbose=False):
    # determination of discrete values
    if discreteExpression==None:
        discreteExpression=discretizeExpressionData(numericalExpression)
    
    # reation of the list of row names of discreteExpression
    if discreteExpression!=None:
        # access to names of the rows of discreteExpression to get a list
        dis_row_names = discreteExpression.index.tolist()

    # determination of GeneList
    if GeneList==None:
        # access to names of the rows of numericalExpression to get a list
        num_row_names = numericalExpression.index.tolist()
        # get genes list witch are not in Tflist
        GeneList = list(set(num_row_names) - set(Tflist))
    
    # determination of maxCoreg
    if maxCoreg==None:
        maxCoreg=len(Tflist)
    
    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION BEFORE STARTING

    '''
        if(  sum(! unique(discreteExpression) %in% -1:1) > 0  ){
            stop("Discrete expression data should only have values in {-1, 0, 1}")}
    '''
    if  any(value not in [-1,0,1] for value in discreteExpression.stack().unique()):
        raise ValueError("Discrete expression data should only have values in {-1, 0, 1}")

    '''    
    if(length(rownames(numericalExpression)) > length(unique(rownames(numericalExpression)))){
        stop("No gene duplicates are allowed in the row.names.")
    }
    ''' 
    if len(num_row_names) > len(set(num_row_names)):
        raise ValueError("No gene duplicates are allowed in the row.names.")

    '''
    if(nrow(numericalExpression) != nrow(discreteExpression) |
    sum(rownames(discreteExpression) != rownames(numericalExpression))>0 ){
        stop("Discrete expression and continuous expression should have the same dimensions and the same rownames (gene/tf names)") }
    '''
    if numericalExpression.shape[0] != discreteExpression.shape[0] | len(num_row_names) > len(dis_row_names):
        raise ValueError("Discrete expression and continuous expression should have the same dimensions and the same rownames (gene/tf names)")

    ''' /!\STOP VERIFICATION : la condition est vraie avec les éléments que j'ai ...'''

    '''
    if(length(intersect(TFlist,rownames(numericalExpression)))<=1 ){
        stop("At least 2 of the provided regulators/transcription factor (TFlist) should be in the rownames in the gene expression matrix")    }
    '''
    if len(set(num_row_names).intersection(Tflist))<=1:
        raise ValueError("At least 2 of the provided regulators/transcription factor (TFlist) should be in the rownames in the gene expression matrix")

    '''
    if(ncol(numericalExpression) > nrow(numericalExpression)){
        warning("Expression data should be in a matrix or data frame with genes in rows and samples in column.")
    }
    '''
    if numericalExpression.shape[1] > numericalExpression.shape[0]:
        raise ValueError("Expression data should be in a matrix or data frame with genes in rows and samples in column.")

    '''
    if(length(intersect(GeneList,rownames(numericalExpression)))==0 ){
        stop("The list of genes (GeneList) should be in the rownames in the gene expression matrix")    }
    '''
    if len(set(num_row_names).intersection(GeneList))==0:
        raise ValueError("The list of genes (GeneList) should be in the rownames in the gene expression matrix")

    '''
    parallel <- match.arg(parallel)
    '''
    parallel_options = ["multicore","no", "snow"]
    if parallel not in parallel_options:
        raise ValueError("The option of parallelism should be \"multicore\", \"no\" or \"snow\"")

    '''
    if(verbose){
        message("Pre-process.")
    '''
    if verbose:
        print("Pre-process.")
    
    
    # select only genes and TF with ones or minus ones in at least minGeneSupport portion of samples
    ''' 
    genesupport = which(apply(abs(discreteExpression), 1 , sum) > (ncol(numericalExpression)*(minGeneSupport)))
    '''
    genes_support = np.where(np.sum(np.abs(discreteExpression), axis=1) > np.shape(numericalExpression)[1] * minGeneSupport)[0]
    
    # to be consistent with R indexing starting at 1
    '''
    discreteExpression=discreteExpression[genesupport,]
    '''
    #genes_support = genes_support + 1  
    discreteExpression=discreteExpression.iloc[genes_support]

    '''
    numericalExpression=numericalExpression[genesupport,]
    '''
    numericalExpression=numericalExpression.iloc[genes_support]

    '''
    TFlist = intersect(rownames(numericalExpression),TFlist)
    '''
    TFlist = list(set(num_row_names).intersection(Tflist))

    '''
    GeneList= intersect(rownames(numericalExpression),GeneList)
    '''
    GeneList= list(set(num_row_names).intersection(GeneList))

    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION AFTER THE SELECTION OF THE GENES AND TF

    ''' /!\STOP VERIFICATION : la condition est vraie avec les éléments que j'ai ...'''
    '''    
    if(length(TFlist)<5){
        stop("Less than 5 of the provided TF are suitable to infer a network. Either provide more TF, more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to select more but less variant TFs.")
    }
    '''
    if len(TFlist)<5:
        raise ValueError("Less than 5 of the provided TF are suitable to infer a network. Either provide more TF, more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to select more but less variant TFs.")

    '''
    if(length(GeneList)==0){
        stop("No genes were suitable to infer regulators. Either provide more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to allow the selection of more but less variant Genes.")
    }
    '''
    if len(GeneList)==0:
        raise ValueError("No genes were suitable to infer regulators. Either provide more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to allow the selection of more but less variant Genes.")
    
    
    # Get all the matrices and datasets needed (gene and tf expression, numerical or discrete)
    # If only one gene is given, R will automatically make a vector. The following make sure this does not happen.
    '''    
    if(length(GeneList)==1){
        geneNumExp= matrix(numericalExpression[GeneList,],nrow=1)
        geneDiscExp= matrix(discreteExpression[GeneList,],nrow=1)
        rownames(geneNumExp)=GeneList
        rownames(geneDiscExp)=GeneList
    }else{
        geneNumExp= numericalExpression[GeneList,]
        geneDiscExp= discreteExpression[GeneList,]
    }
    '''
    if len(GeneList)==1:
        geneNumExp = numericalExpression.loc[GeneList]
        geneDiscExp= discreteExpression.loc[GeneList]
        geneNumExp.index=GeneList
        geneDiscExp.index=GeneList
    else: # le else a-t-il un intérêt
        geneNumExp= numericalExpression.loc[GeneList]
        geneDiscExp= discreteExpression.loc[GeneList]
    
    '''
    regNumExp= numericalExpression[TFlist,]
    regDiscExp= discreteExpression[TFlist,]
    '''
    regNumExp= numericalExpression.loc[TFlist]
    regDiscExp= discreteExpression.loc[TFlist]
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## TRANSFORMING ALL DISCRETE DATA INTO TRANSACTIONs
    # To run apriori, the discrete data must be binary. So, the discrete data is simply becoming two concatenated binary matrix
    # first n samples are positive expression values, then all negative values.
    '''
    posSamples = 1:ncol(discreteExpression)
    negSamples= (ncol(discreteExpression) +1):(ncol(discreteExpression) *2)
    regBitData =cbind(regDiscExp==+1 , regDiscExp== -1)
    transRegBitData= as(t(regBitData),"transactions")
    '''
    '''                                                 /!\ PAUSE DANS LES VERIFICATIONS LE 12 04 2024 '''
    '''
    if(verbose){
        message("Mining coregulator ...")
    }
    '''
    posSamples = range(1,discreteExpression.shape[1]) # non utilisé

    negSamples= range(discreteExpression.shape[1]+1,discreteExpression.shape[1] *2) # non utilisé

    regBitData =np.column_stack((regDiscExp == +1, regDiscExp == -1))

    transRegBitData= regBitData.T

    if verbose:
        print("Mining coregulator ...")
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## MINING FOR FREQUENT COREGULATORS
    # using apriori instead of eclat. testing may be required for possible speed improvement.
    '''
    miningFunction=apriori
    transitemfreq =suppressWarnings(miningFunction(transRegBitData,parameter=list(support = minGeneSupport/2,maxlen=1,target="frequent itemsets")
    ,control=list(verbose=FALSE)))
    if(maxCoreg > 1){
        transitemfreq=c(transitemfreq,suppressWarnings(miningFunction(transRegBitData,parameter=list(support =minCoregSupport/2,minlen=2,maxlen=maxCoreg,target="closed frequent itemsets")
        ,control=list(verbose=FALSE))))
    }
    coregs =as(slot(transitemfreq,"items"),"list")
    '''
    te = TransactionEncoder()
    te_ary = te.fit(transRegBitData).transform(transRegBitData)
    transRegBitData_df = pd.DataFrame(te_ary, columns=te.columns_)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        transitemfreq=apriori(transRegBitData_df, min_support=0.5, use_colnames=True, max_len=1, verbose=0, low_memory=False)
        # low_memory=True should only be used for large dataset if memory resources are limited 
    
        if maxCoreg > 1:
            result=apriori(transRegBitData_df, min_support=minCoregSupport/2, use_colnames=True, max_len=maxCoreg, verbose=0, low_memory=False)
            transitemfreq.extend(result)

        coregs = transitemfreq['itemsets'].tolist()

    '''             /!\ QUESTION ????
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
              +" target genes, "+str(len(TFlist))
              +" regulators and a total of coregulator sets "+str(len(coregs))
              +" sets of potential co-regulators.\nSearch parameters :\n"+"Maximum size of co-regulator sets : "+str(maxCoreg)
              +"\nNumber of putative GRN per gene : "+str(nGRN)+"\nMinimum number of differentially expressed samples to select a single gene : "+str(minGeneSupport)+
              "\nMinimum number of differentially expressed samples to select a set of co-regulator : "+str(minCoregSupport))
        
        print("Mining GRN ...")
    
    '''
    result=data.frame()
    gotNet=FALSE
    '''
    result=np.dataframe()
    gotNet=False

    #just because it's easier toadd here 5% and remove it at the first line in the while loop, where it needs to be decrementale in case no GRNs are found
    '''searchThresh=  1/((1/searchThresh)-1)'''
    searchThresh=  1/((1/searchThresh)-1)
    
    # In very large datasets of very heterogeneous samples (such as the large collection of unrelated cell lines ...)
    # It is possible that no GRN can be fitted with stringent threshold (usually 50%) and that no GRN is found.
    # In case this happens, the threshold is decremented step by step and if no network is found at 10%, then none can be found ...
    '''
    while(searchThresh >= 0.05 & !gotNet )
    {
        searchThresh =  1/((1/searchThresh)+1)
        if(parallel =="multicore" & length(GeneList)>1 & getOption("mc.cores", 2L) > 1)
        {
            result =mclapply(GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else if(parallel =="snow" & !is.null(cluster) & length(GeneList)>1){
            result =parLapply(cluster,GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else if( length(GeneList)>1){
            result =lapply(GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else{
            result=oneGeneHLICORN(GeneList,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh ,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }
    }'''
    while searchThresh >= 0.05 & gotNet==True :
        # decrements the search threshold in case nothing is found
        #(can be the case for VERY large datasets for which it can be hard to find regulators with 50% of matching +1 and -1)
        searchThresh =  1/((1/searchThresh)+1)
        #running hlicorn for each gene in a multithread way if needed.
        
        # get available cares on the machine and define how much to use
        available_cores_nb = multiprocessing.cpu_count()
        using_processus = max(1, available_cores_nb - 1)

        def process_gene(gene):
            return oneGeneHLICORN(gene, geneDiscExp, regDiscExp, coregs, transitemfreq, transRegBitData, searchThresh, genexp=geneNumExp, regnexp=regNumExp, nresult=nGRN)

        if parallel =="multicore" & len(GeneList)>1 & using_processus > 1:
            '''             /!\ QUESTION ????
            ProcessPoolExecutor est un composant de la bibliothèque standard de Python qui permet d'exécuter des fonctions de manière asynchrone dans des processus parallèles.
            Il crée un pool de processus dans lequel les fonctions peuvent être exécutées de manière concurrente.

            L'avantage de ProcessPoolExecutor par rapport à ThreadPoolExecutor est qu'il utilise des processus plutôt que des threads.
            Cela signifie que chaque fonction est exécutée dans son propre processus, ce qui peut être bénéfique pour les tâches qui impliquent des opérations CPU-bound,
            car elles peuvent être parallélisées sur plusieurs cœurs de CPU.
            '''
            
            '''
            Version avec joblib
            
            import multiprocessing
            
            processes = [] # facultatif, uniquement pour récupérer la liste des processus qui se sont exécutés
            for _ in GeneList:
                p = multiprocessing.Process(target=process_gene,args=[GeneList])
                p.start()
                processes.append(p)
            
            for process in processes:
                process.join() # pour forcer le programme à attendre que tous les process soient fini avec de continuer
                
            peu importe le nombre de coeurs disponibles, les taches seront distribuées en fonction des ressources disponibles.
            '''
            
            # create a pool with the number of processus wanted
            with ProcessPoolExecutor() as executor:
                # the results are in the list "results"
                results = list(executor.map(process_gene, GeneList))
                # je ne suis pas sûre de devoir laisser list() quand on utilise map(). map() garde une liste en mémoire contrairement à imap et imap_unordered
                
            '''
            variante : results = [executor.submit(process_gene(gene)) for gene in GeneList]
            à privilégier par rapport à multiprocessing ?!
            mais il vaut mieux utiliser map qu'une boucle pour que les process s'exécutent dans l'ordre de la liste
            '''
            
            
            '''             /!\ QUESTION ????
            version avec les Threads :

            from concurrent.futures import ThreadPoolExecutor

            with ThreadPoolExecutor() as executor:
                results = list(executor.map(process_gene, GeneList))
            
            Il faut tester les méthodes ProcessPoolExecutor et ThreadPoolExecutor car dans certains cas, l'une est plus rapide que l'autre
            '''
            
            '''
            Version avec joblib :
                
            from joblib import Parallel, delayed, parallel_config

            with parallel_config(backend="loky", inner_max_num_threads=2):
                results = Parallel(n_jobs=using_processus)(delayed(process_gene)(gene) for gene in GeneList)
            '''
            gotNet=True

        elif parallel =="snow" & cluster!=None & len(GeneList)>1:
            with multiprocessing.Pool() as pool:
                results = [pool.apply_async(process_gene, (gene,)) for gene in GeneList]
                results = [res.get() for res in results]
            '''             /!\ QUESTION ????
            Pool() crée un pool de processus avec un nombre de processus par défaut (nombre de cœurs du CPU).
            pool.apply_async est utilisé pour exécuter de manière asynchrone la fonction process_gene pour chaque élément de GeneList.
            res.get() est utilisé pour obtenir les résultats de chaque tâche asynchrone.
            
                        /!\ MULTIPROCESSING.POOL
            à utiliser avec une protection du main
            if __name__ == '__main__':
            '''
            gotNet=True

        elif len(GeneList)>1:
            # comprehension of the list to apply prosses_gene to each gene of the list
            results = [process_gene(gene) for gene in GeneList]
                        
            gotNet=True

        else:
            result=oneGeneHLICORN(GeneList,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)

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
        if results.shape[1] >= 3 & results.shape[0] >0:
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
    if results.shape[0]==0 | results.shape[1] <3:
        raise ValueError("Something went wrong. No GRN found.")

    results.iloc[:, 0] = None

    '''         /!\ fin à éclaircir
        sigrns = coregnet(result)
    sigrns@inferenceParameters=list(minGeneSupport=minGeneSupport,maxCoreg=maxCoreg,minCoregSupport = minCoregSupport,searchThresh=searchThresh,nGRN=nGRN)
        return(sigrns)'''
        