import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from functools import reduce

'''
discretizeExpressionData = function(numericalExpression,threshold=NULL,refSamples=NULL,standardDeviationThreshold=1)
'''
def discretizeExpressionData(numericalExpression,threshold=None,refSamples=None,standardDeviationThreshold=1):
    '''  
    numericalExpression=as.matrix(numericalExpression) --> pas besoin de le transcrire en python
    
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

        centered.index = numericalExpression.index

        #colnames(centered) = setdiff(colnames(numericalExpression),refSamples)
        columns_to_exclude = set(refSamples)
        columns_to_use = reduce(lambda left, right: left.difference(right), [numericalExpression.columns, columns_to_exclude])
        centered.columns = columns_to_use

    elif np.min(numericalExpression) >= 0: #  means that it's raw (log or not) data
        #centered  =t(scale(t(numericalExpression),scale=FALSE))
        centered = numericalExpression.sub(numericalExpression.mean(axis=1), axis=0)
        centered.index = numericalExpression.index
        centered.columns = numericalExpression.columns
        
    else:
        centered[pd.isna(numericalExpression)] = 0

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