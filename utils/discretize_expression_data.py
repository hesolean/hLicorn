# librabies
import pandas as pd
from sklearn.preprocessing import StandardScaler


def discretize_expression_data(numerical_expression, refSamples=None) :
    #########################################
    # take off all the NaN values
    numerical_expression[pd.isna(numerical_expression)]=0

    # create object before scaling
    scaler=StandardScaler(with_mean=True, with_std=False)

    # if refSample exist, scaling according to refSample values else scaling according to every values
    if refSamples != None :
        scaler.fit(numerical_expression[refSamples])
    else :
        scaler.fit(numerical_expression)
    scaled_data=scaler.transform(numerical_expression)
    
    # discretisize the values
    scaled_data[scaled_data <= -1]=-1
    scaled_data[scaled_data >= 1]=1
    scaled_data[(scaled_data > -1) & (scaled_data < 1)]=0

    # transform in integer values
    scaled_data=scaled_data.astype(int)

    # take off all the NaN values
    scaled_data[pd.isna(scaled_data)]=0

    # transform a numpy.ndArray in pandas dataframe
    scaled_data=pd.DataFrame(scaled_data)
    
    # rename rows and columns
    scaled_data.index=numerical_expression.index
    scaled_data.columns=numerical_expression.columns

    # return dataframe pandas
    return(scaled_data)
    