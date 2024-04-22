# librabies
import pandas as pd

# secondary functions
from utils.hLicorn import hLICORN

def main():
    numericalExpression = pd.read_csv('CIT.csv')

    # rownames as index to exit str in the dataframe
    numericalExpression.set_index(numericalExpression.columns[0], inplace=True)

    Tflist = pd.read_csv('HumanTF.csv')

    try:
       hLICORN(numericalExpression,Tflist)
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()