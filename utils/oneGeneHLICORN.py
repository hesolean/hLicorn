import numpy as np
import ctypes

from utils.eand import eand


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
    /!\ JE N'AI PAS TROUVE A QUOI CORRESPOND 'support'
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
    coact[coactnames.iloc[g, :]==("")]=""

    corep_sorted = corep.apply(sorted, axis=1)
    corep_pasted = corep_sorted.apply(lambda x: ' '.join(map(str, x)))
    corepnames = corep_pasted.unique()

    corep=corepnames.split()
    corep[corepnames.iloc[g, :]==("")]=""

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

    '''on appelle le bout de programme en C, A VERIFIER
    x= .C("combnLicorn",
    as.integer(t(coactexp)),as.integer(length(coact)),#expression and number of coactivators
    as.integer(t(corepexp)),as.integer(length(corep)),#expression and number of corepressors
    as.integer(geneDiscExp[g,]),      as.integer(ncol(geneDiscExp)),    #expression of gene, number of samples
    as.double(rep(-1,(length(coact))*(length(corep)))),# vector to store MAE results
    as.integer(rep(-1,(length(coact))*(length(corep)))),# vector to store index of coactivator
    as.integer(rep(-1,(length(coact))*(length(corep))))# vector to store index of corepressor
    )
    '''
    

    # Chargement la bibliothèque partagée contenant la fonction C, /!\ AJUSTER LE CHEMIN
    lib = ctypes.CDLL("./library.so")

    # Définition des types des arguments de la fonction C
    lib.combnLicorn.argtypes = [
        ctypes.POINTER(ctypes.c_int),  # coactexp
        ctypes.POINTER(ctypes.c_int),  # ncoacts
        ctypes.POINTER(ctypes.c_int),  # corepexp
        ctypes.POINTER(ctypes.c_int),  # ncoreps
        ctypes.POINTER(ctypes.c_int),  # gexp
        ctypes.POINTER(ctypes.c_int),  # nsamples
        ctypes.POINTER(ctypes.c_double),  # result
        ctypes.POINTER(ctypes.c_int),  # iact
        ctypes.POINTER(ctypes.c_int),  # irep
    ]

    # Définition des valeurs des arguments
    coactexp = np.array(coactexp)
    ncoacts = np.array([len(coactexp)], dtype=np.int32)
    corepexp = np.array(corepexp)
    ncoreps = np.array([len(corepexp)], dtype=np.int32)
    gexp = np.array(geneDiscExp)  # A VERIFIER
    nsamples = np.array([len(g)], dtype=np.int32)
    result = np.zeros((ncoacts[0], ncoreps[0]), dtype=np.float64)
    iact = np.zeros((ncoacts[0], ncoreps[0]), dtype=np.int32)
    irep = np.zeros((ncoacts[0], ncoreps[0]), dtype=np.int32)

    # Appel de la fonction C avec les arguments appropriés
    lib.combnLicorn(
        coactexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoacts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        corepexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoreps.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        gexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        nsamples.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        iact.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        irep.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )

    # result, iact, irep contiennent les résultats retournés par la fonction C

    
    # bad index will store all bad comparisons (could be done before computing .. right?)
    # then, no intersection between act and rep (includes both empty

    ''' A PARTIR DE LA, JE N'AI PAS COMPRIS EXACTEMENT L'INTERET DONC A REVOIR 
    // on combine les colonnes 8 et 9 qui viennent de x avec la longueur de l'intersection entre un coactivateur et un coinhibiteur
    pour finir on renvoie les indices des valeurs égales à 0
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
