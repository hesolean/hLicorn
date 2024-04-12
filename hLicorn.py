import pandas as pd
import numpy as np
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori
import warnings
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from sklearn.preprocessing import StandardScaler
from functools import reduce

'''
hLICORN=function( numericalExpression,discreteExpression=discretizeExpressionData(numericalExpression)
, TFlist, GeneList=setdiff(rownames(numericalExpression),TFlist),parallel = c("multicore","no", "snow"),cluster=NULL,
minGeneSupport=0.1,minCoregSupport = 0.1,maxCoreg=length(TFlist),
searchThresh=1/3,nGRN=100,verbose=FALSE)
'''    
def hLICORN(numericalExpression,Tflist,discreteExpression=None,GeneList=None,parallel="no",cluster=None,minGeneSupport=0.1,minCoregSupport = 0.1,maxCoreg=None,searchThresh=1/3,nGRN=100,verbose=False):
    # determination of discrete values
    if discreteExpression==None:
        discreteExpression=discretizeExpressionData(numericalExpression)
    
    # reation of the list of row names of discreteExpression
    if discreteExpression!=None:
        # tranformation of discreteExpression in df
        dis_exp_df = pd.DataFrame(discreteExpression)
        
        # access to names of the rows of discreteExpression to get a list
        dis_row_names = dis_exp_df.iloc[:, 0].tolist()

    # determination of GeneList
    if GeneList==None:
        # access to names of the rows of numericalExpression to get a list
        num_row_names = numericalExpression.iloc[:, 0].tolist()
        # get genes list witch are not in Tflist
        GeneList=num_row_names.difference(Tflist)
    
    # determination of maxCoreg
    if maxCoreg==None:
        maxCoreg=len(TFlist)
    
    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION BEFORE STARTING

    '''
        if(  sum(! unique(discreteExpression) %in% -1:1) > 0  ){
            stop("Discrete expression data should only have values in {-1, 0, 1}")}
    '''
    if  any(value not in [-1,0,1] for value in discreteExpression.unique()):
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
    indices_genes_support = np.where(np.sum(np.abs(discreteExpression), axis=1) > np.shape(numericalExpression)[1] * minGeneSupport)[0]
    
    # to be consistent with R indexing starting at 1
    '''
    discreteExpression=discreteExpression[genesupport,]
    '''
    genes_support = indices_genes_support + 1  
    discreteExpression=discreteExpression[genes_support,:]

    '''
    numericalExpression=numericalExpression[genesupport,]
    '''
    numericalExpression=numericalExpression[genes_support,:]

    '''
    TFlist = intersect(rownames(numericalExpression),TFlist)
    '''
    TFlist = set(num_row_names).intersection(Tflist)

    '''
    GeneList= intersect(rownames(numericalExpression),GeneList)
    '''
    GeneList= set(num_row_names).intersection(GeneList)

    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION AFTER THE SELECTION OF THE GENES AND TF
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
        geneNumExp = np.matrix(numericalExpression[GeneList, :]).T
        geneDiscExp= np.matrix(discreteExpression[GeneList, :]).T
        geneNumExp.iloc[:, 0]=GeneList
        geneDiscExp.iloc[:, 0]=GeneList
    else:
        geneNumExp= numericalExpression[GeneList, :]
        geneDiscExp= discreteExpression[GeneList, :]
    
    '''
    regNumExp= numericalExpression[TFlist,]
    regDiscExp= discreteExpression[TFlist,]
    '''
    regNumExp= numericalExpression[TFlist, :]
    regDiscExp= discreteExpression[TFlist, :]
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## TRANSFORMING ALL DISCRETE DATA INTO TRANSACTIONs
    # To run apriori, the discrete data must be binary. So, the discrete data is simply becoming two concatenated binary matrix
    # first n samples are positive expression values, then all negative values.
    '''
    posSamples = 1:ncol(discreteExpression)
    negSamples= (ncol(discreteExpression) +1):(ncol(discreteExpression) *2)
    regBitData =cbind(regDiscExp==+1 , regDiscExp== -1)
    transRegBitData= as(t(regBitData),"transactions")
    
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
        print("Learning a Co-Regulatory network for:\n"+str(len(GeneList))+" target genes, "+str(len(TFlist))+" regulators and a total of coregulator sets "+str(len(coregs))+" sets of potential co-regulators.\nSearch parameters :\n"+"Maximum size of co-regulator sets : "+str(maxCoreg)+"\nNumber of putative GRN per gene : "+str(nGRN)+"\nMinimum number of differentially expressed samples to select a single gene : "+str(minGeneSupport)+"\nMinimum number of differentially expressed samples to select a set of co-regulator : "+str(minCoregSupport))
        
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
            
            # create a pool with the number of processus wanted
            with ProcessPoolExecutor() as executor:
                # the results are in the list "results"
                results = list(executor.map(process_gene, GeneList))
            '''             /!\ QUESTION ????
            version avec les Threads :

            from concurrent.futures import ThreadPoolExecutor

            with ThreadPoolExecutor() as executor:
                results = list(executor.map(process_gene, GeneList))
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
            gotNet=FALSE
        
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
        
'''fin hLICORN'''



'''
discretizeExpressionData = function(numericalExpression,threshold=NULL,refSamples=NULL,standardDeviationThreshold=1)
'''
def discretizeExpressionData(numericalExpression,threshold=None,refSamples=None,standardDeviationThreshold=1):
    '''  
    numericalExpression=as.matrix(numericalExpression)
    '''
    numericalExpression = numericalExpression.values  
  
    '''
    if(!is.null(refSamples) ){
        refmeans = apply(numericalExpression[,refSamples],1,mean)
        centered  =t(scale(t(numericalExpression[,setdiff(colnames(numericalExpression),refSamples)]),scale=FALSE,center=refmeans))  
        rownames(centered) = rownames(numericalExpression)
        colnames(centered) = setdiff(colnames(numericalExpression),refSamples)
    }else if(min(matrix(numericalExpression) ) >= 0){#  means that it's raw (log or not) data
        centered  =t(scale(t(numericalExpression),scale=FALSE))  
        dimnames(centered) = dimnames(numericalExpression) 
    }else{
        centered=numericalExpression
        centered[which(is.nan(centered))]=0
    }
    '''
    if refSamples!=None:
        refmeans = numericalExpression[refSamples].apply(np.mean, axis=1)

        #centered =t(scale(t(numericalExpression[,setdiff(colnames(numericalExpression),refSamples)]),scale=FALSE,center=refmeans))  
        # columns to center and scaling
        selected_columns = numericalExpression.columns.difference(refSamples)

        # mean of the columns
        refmeans = numericalExpression[refSamples].mean(axis=1)

        # transpose, scalling and retranspose
        transpose_data = numericalExpression[selected_columns].T
        scaler = StandardScaler(with_mean=False, with_std=False)
        scaled_data = scaler.fit_transform(transpose_data)
        retranspose_data = pd.DataFrame(scaled_data.T, columns=selected_columns)

        # center the data
        centered_data = retranspose_data.sub(refmeans, axis=0)

        # finale transposition
        centered = centered_data.T

        centered.iloc[:, 0] = numericalExpression.iloc[:, 0]

        #colnames(centered) = setdiff(colnames(numericalExpression),refSamples)
        columns_to_exclude = set(refSamples)
        columns_to_use = reduce(lambda left, right: left.difference(right), [numericalExpression.columns, columns_to_exclude])
        centered.columns = columns_to_use

    elif min(numericalExpression) >= 0: #  means that it's raw (log or not) data
        #centered  =t(scale(t(numericalExpression),scale=FALSE))
        centered = numericalExpression.sub(numericalExpression.mean(axis=1), axis=0).T
        centered.iloc[:, 0] = numericalExpression.iloc[:, 0]
        centered.columns = numericalExpression.columns
        
    else:
        centered=numericalExpression.fillna(0, inplace=True)

    '''
    if(is.null(threshold)){
        threshold=sd(centered)*standardDeviationThreshold
    }
    nco = ncol(centered)
    nro = nrow(centered)
    discreteExpression =matrix(  as.integer( centered >= threshold  ) + (-as.integer(centered <= (- threshold) )),nrow=nro,ncol=nco)
    dimnames(discreteExpression) = dimnames(centered)
    return(discreteExpression)
    }
    '''
    if threshold!=None:
        centered_std = centered.std()
        threshold = centered_std * standardDeviationThreshold

    '''         /!\ A VOIR
    transcrit en python mais à priori pas besoin dans le calcul de discreteExpression
    nco = centered.shape[1]
    nro = centered.shape[0]
    '''
    #discreteExpression =matrix(  as.integer( centered >= threshold  ) + (-as.integer(centered <= (- threshold) )),nrow=nro,ncol=nco)
    above_threshold = centered >= threshold
    below_neg_threshold = centered <= -threshold

    above_threshold_int = above_threshold.astype(int)
    below_neg_threshold_int = below_neg_threshold.astype(int)

    discreteExpression = above_threshold_int + (-below_neg_threshold_int)
    discreteExpression.iloc[:, 0] = centered.iloc[:, 0]
    discreteExpression.column = centered.column

    return(discreteExpression)

    '''Fin discretizeExpressionData'''


'''oneGeneHLICORN = function(g,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,regnexp,genexp,nresult)'''
def oneGeneHLICORN(g,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,regnexp,genexp,nresult):
    '''{
    shift=ncol(geneDiscExp)
    '''
    shift=geneDiscExp.shape[1]

    '''
    pos =which( geneDiscExp[g,]==1)
    neg=which( geneDiscExp[g,]== -1)+shift
    '''
    #sample index with the target gene at 1 or -1
    pos = np.where(geneDiscExp.iloc[g, :] == 1)[0]
    # negative samples are shifted because we are using a binary matrix with true false for ones in the first part
    # and true false for -1 in the second part
    neg = np.where(geneDiscExp.iloc[g, :] == -1)[0]+shift
    
    '''   
    coact=coregs[which(support(transitemfreq, transRegBitData[c(pos,neg)])>= searchThresh )]
    pos =pos+shift
    neg=neg - shift
    corep=coregs[which(support(transitemfreq, transRegBitData[c(pos,neg)]) >= searchThresh )]
    '''
    # select all the coregulators with a support of 50% minimum only in the samples with the target gene at ones or minus ones
    # indices for which the threshold is reached, then we select the elements
    coact_indices = np.where(support(transitemfreq, transRegBitData.iloc[:, pos + neg]) >= searchThresh)[0]
    coact = coregs.iloc[coact_indices]

    pos=pos+shift
    neg=neg-shift
    corep_indices = np.where(support(transitemfreq, transRegBitData.iloc[:, pos + neg]) >= searchThresh)[0]
    corep = coregs.iloc[corep_indices]
    
    '''
    corep = c(corep,list(""))
    coact = c(coact,list(""))
    '''
    # add empty coregulators to have the possibility to only have ativators or inhibitors
    corep.append("")
    coact.append("")
    
    '''    
    coactnames =unique(sapply(lapply(coact,sort),paste,collapse=" "))
    coact=strsplit(coactnames," ")
    coact[[which(coactnames=="")]]=""
    corepnames =unique(sapply(lapply(corep,sort),paste,collapse=" "))
    corep=strsplit(corepnames," ")
    corep[[which(corepnames=="")]]=""
    '''
    # to have unique coregulators and a single vector of coreg (not a list)
    #coactnames =unique(sapply(lapply(coact,sort),paste,collapse=" "))
    coact_sorted = coact.apply(sorted, axis=1)
    coact_pasted = coact_sorted.apply(lambda x: ' '.join(map(str, x)))
    coactnames = coact_pasted.unique()

    coact=coactnames.split()
    coact[coactnames.index("")]=""

    corep_sorted = corep.apply(sorted, axis=1)
    corep_pasted = corep_sorted.apply(lambda x: ' '.join(map(str, x)))
    corepnames = corep_pasted.unique()

    corep=corepnames.split()
    corep[corepnames.index("")]=""

    '''   
    coactexp = eand(coact,regDiscExp)
    corepexp = as.integer(eand(corep,regDiscExp))
    '''
    # merge merge expression of coregulator and corepressor
    coactexp = eand(coact,regDiscExp)
    corepexp = int(eand(corep,regDiscExp))

    '''
    corepexp[which(corepexp ==1)] = 2
    '''
    # active inhibitor has a stronger impact than naything else in licorn:
    corepexp.replace(1, 2, inplace=True)    

    '''on appelle le bout de programme en C, mis de coté pour l'instant...
    x= .C("combnLicorn",
    as.integer(t(coactexp)),as.integer(length(coact)),#expression and number of coactivators
    as.integer(t(corepexp)),as.integer(length(corep)),#expression and number of corepressors
    as.integer(geneDiscExp[g,]),      as.integer(ncol(geneDiscExp)),    #expression of gene, number of samples
    as.double(rep(-1,(length(coact))*(length(corep)))),# vector to store MAE results
    as.integer(rep(-1,(length(coact))*(length(corep)))),# vector to store index of coactivator
    as.integer(rep(-1,(length(coact))*(length(corep))))# vector to store index of corepressor
    )
    '''
    # bad index will store all bad comparisons (could be done before computing .. right?)
    # then, no intersection between act and rep (includes both empty

    ''' A PARTIR DE LA, JE N'AI PAS COMPRIS EXACTEMENT L'INTERET DONC A REVOIR 
    // on combine les colonnes 8 et 9 qui viennent de x avec la longueur de l'intersection entre un coactivateur et un coinhibiteur
    pour finir on renvie les indices des valeurs égales à 0
    goodindex=which(apply(cbind(x[[8]],x[[9]]),1,function(y){
        return(length(intersect( coact[[y[1]]],corep[[y[2]]] )))
    })==0)
    
    
    selact = coactnames[x[[8]][goodindex]]
    selrep = corepnames[x[[9]][goodindex]]
    
    # all empty set of coregulators are set to NA
    selact[which(selact=="")]=NA
    selrep[which(selrep=="")]=NA
    // calcul de performence du modèle de regression Mean Absolute Error
    mae = x[[7]][goodindex]
    if(!is.na(nresult)){
        # get 100 first ranks, if ties, might get more ...
        bestindex= which(rank(mae,ties.method="min")<=nresult)
        GRN = data.frame("Target"=rep(g, length(bestindex)),"coact"=selact[bestindex],   "corep"=selrep[bestindex] ,stringsAsFactors=FALSE)      
    }else{
        bestindex= which(rank(mae,ties.method="min")==1)
        GRN = data.frame("Target"=rep(g, length(bestindex)),"coact"=selact[bestindex],   "corep"=selrep[bestindex] ,stringsAsFactors=FALSE)      
        
    }
  
    # if no grn are found return NULL
    if(nrow(GRN)==0){
        return(NULL)
    }
    
    if(is.matrix(GRN) | is.data.frame(GRN)){
        linearmodels=.getEntry(apply(GRN,1,.fitGRN,genexp=t(genexp),regexp=t(regnexp),permut=FALSE),"numscores")
    }else{
        #    linearmodels=.linearCoregulationBootstrap(as.character(GRN),genexp=gexp,regnexp=regnexp,numBootstrap=numBootstrap)
        linearmodels=.fitGRN(as.character(GRN),genexp=t(genexp),regexp=t(regnexp),permut=FALSE)$numscores
    }
    
    numscores=data.frame(t(linearmodels),stringsAsFactors = FALSE)
    
    numscores[,3]=as.numeric(numscores[,3])
    numscores[,4]=as.numeric(numscores[,4])
    numscores[,5]=as.numeric(numscores[,5])
    numscores[,6]=as.numeric(numscores[,6])
    colnames(GRN)=c("Target","coact","corep")
    
    return(data.frame(GRN,numscores,stringsAsFactors = FALSE))
    
    }'''
    '''Fin oneGeneHLICORN'''

'''eand = function(coact,regDiscExp,multip=1)'''
def eand(coact,regDiscExp,multip=1):
    '''
    do.call(rbind,lapply(coact,function(co){
        if(co[1] == ""){
            return(rep(0,ncol(regDiscExp)))
        }else if(length(co) ==1){return(regDiscExp[co,])}
        n =length(co)
        x=apply(regDiscExp[co,],2,sum)
        y=x
        x[1:length(x)]=0
        x[which(y== - n )]=- multip
        x[which(y == n)] = multip
        return(x)
    }))
    '''
    def function_co(co):
        if co[1] == "":
            return np.zeros(regDiscExp.shape[1])
        elif len(co)==1:
            return regDiscExp[co, :]
        n=len(co)
        x = regDiscExp.iloc[co, :].sum(axis=0)
        y=x
        x[1:len(x)]=0
        x.loc[y[y==-n].index]=-multip
        x.loc[y[y==n].index]=multip
        return x

    concatened_df = pd.concat(map(lambda co: function_co(co),coact))
    return concatened_df

''' Fin eand '''

'''         MAIN            '''
numericalExpression = pd.read_csv('CIT.csv')
try:
    hLICORN()
except ValueError as e:
    print(e)