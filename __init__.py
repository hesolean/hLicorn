# librabies
import pandas as pd

# secondary functions
from utils.hLicorn import hLICORN

def main():
    numerical_expression = pd.read_csv('CIT.csv')

    # rownames as index to exit str in the dataframe
    numerical_expression.set_index(numerical_expression.columns[0], inplace=True)

    tf_list = pd.read_csv('HumanTF.csv')

    try:
       hLICORN(numerical_expression,tf_list, parallel="no", max_coreg=3)
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()