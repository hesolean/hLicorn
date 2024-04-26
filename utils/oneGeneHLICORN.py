# librabies
import numpy as np
from mlxtend.frequent_patterns import association_rules

# secondary functions
from utils.eand import eand
from utils.C_call import C_call


'''oneGeneHLICORN = function(g,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,regnexp,genexp,nresult)'''
def oneGeneHLICORN(g,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,regnexp,genexp,nresult):
    '''{
    shift=ncol(geneDiscExp)
    '''
    shift=geneDiscExp.shape[1]
    print(f"shift : {shift}")

    '''
    pos =which( geneDiscExp[g,]==1)
    neg=which( geneDiscExp[g,]== -1)+shift
    '''
    #sample index with the target gene at 1 or -1
    pos = np.where(geneDiscExp.loc[g, :] == 1)[0]
    # negative samples are shifted because we are using a binary matrix with true false for ones in the first part transRegBitData(pos,neg)
    # and true false for -1 in the second part
    neg = np.where(geneDiscExp.loc[g, :] == -1)[0]+shift
    print(f"pos de taille : ", len(pos)," neg de taille : ", len(neg))
    print(f"pos : {pos}")
    print(f"neg : {neg}")
    '''
    coact=coregs[which(support(transitemfreq, transRegBitData[c(pos,neg)])>= searchThresh )]
    pos =pos+shift
    neg=neg - shift
    corep=coregs[which(support(transitemfreq, transRegBitData[c(pos,neg)]) >= searchThresh )] # on relance le calcul d'itemset frequent mais sur une liste restreinte car que pour le gène en cours d'étude
    '''
    # select all the coregulators with a support of 50% minimum only in the samples with the target gene at ones or minus ones
    # indices for which the threshold is reached, then we select the elements
    coact = association_rules(transitemfreq, metric="support", min_threshold=searchThresh) # pas besoin et support est une colonne du df donc dans le calcul des itemse, on peut filtrer directement les support avec le searchThresh et on le fait en même temps que le calcul des itemset
    # !!!!!!!!!!!!!! aller sur le doc issue du github de mlxtend
    print(f"coact : {coact}")
    pos=pos+shift # on garde le fait de tester les pos neg puis les neg pos
    neg=neg-shift
    #corep =   liste des itemset qui sortent du df de frequent_itemset, il faut le jeu de données liée au gène étudié, transRegBitData(pos,neg)
    print(f"corep : {corep}")

    '''
    corep = c(corep,list(""))
    coact = c(coact,list(""))
    '''
    # add empty coregulators to have the possibility to only have ativators or inhibitors
    corep.loc[len(corep)] = ""
    coact.loc[len(coact)] = ""

    # to have unique coregulators and a single vector of coreg (not a list)
    '''   
    juste les listes des noms triées 
    coactnames =unique(sapply(lapply(coact,sort),paste,collapse=" "))
    coact=strsplit(coactnames," ")
    coact[[which(coactnames=="")]]=""
    corepnames =unique(sapply(lapply(corep,sort),paste,collapse=" "))
    corep=strsplit(corepnames," ")
    corep[[which(corepnames=="")]]=""
    '''
    coact = coact.drop_duplicates() # vérifier si utile
    corep = corep.drop_duplicates()

    # merge expression of coregulator and corepressor
    print(f"regDiscExp : {regDiscExp}")

    coactexp = eand(coact, regDiscExp)
    corepexp = eand(corep, regDiscExp)

    print(f"coactexp : {coactexp}")
    print(f"corepexp : {corepexp}")

    # active inhibitor has a stronger impact than naything else in licorn:
    '''
    corepexp[which(corepexp ==1)] = 2
    '''
    #C_call(coactexp,corepexp,geneDiscExp)

    '''
    # bad index will store all bad comparisons (could be done before computing .. right?)
    # then, no intersection between act and rep (includes both empty
    goodindex=which(apply(cbind(x[[8]],x[[9]]),1,function(y){
        return(length(intersect( coact[[y[1]]],corep[[y[2]]] )))
    })==0)
    
    
    selact = coactnames[x[[8]][goodindex]]
    selrep = corepnames[x[[9]][goodindex]]
    
    # all emty set of coregulators are set to NA
    selact[which(selact=="")]=NA
    selrep[which(selrep=="")]=NA
    
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
