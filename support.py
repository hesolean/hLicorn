def support(items, test_df):
	freq_items = pd.DataFrame()
	supports = []
	for itemset in items.itemsets:
		freq = test_df.shape[0]
		for item in itemset:
			sup = 0
			if item in test_df.columns:	sup = test_df[test_df[item] == True][item].count()
			if sup<freq :	freq=sup
		print(str(itemset), " : " , str(freq))
		supports.append(freq/test_df.shape[0])

	freq_items["support"]=supports
	freq_items["items"]=items.itemsets
	
	return freq_items




############# juste pour tester la fonction
dataset = [['Milk', 'Onion', 'Nutmeg', 'Kidney Beans', 'Eggs', 'Yogurt'],
           ['Dill', 'Onion', 'Nutmeg', 'Kidney Beans', 'Eggs', 'Yogurt'],
           ['Milk', 'Apple', 'Kidney Beans', 'Eggs'],
           ['Milk', 'Unicorn', 'Corn', 'Kidney Beans', 'Yogurt'],
           ['Dill', 'Onion', 'Nutmeg', 'Kidney Beans', 'Eggs', 'Yogurt'],
           ['Milk', 'Apple', 'Kidney Beans', 'Eggs'],
           ['Milk', 'Unicorn', 'Corn', 'Kidney Beans', 'Yogurt'],
           ['Corn', 'Onion', 'Onion', 'Kidney Beans', 'Ice cream', 'Eggs']]
		   
		   

import pandas as pd
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori


te = TransactionEncoder()
te_ary = te.fit(dataset).transform(dataset)
df = pd.DataFrame(te_ary, columns=te.columns_)
print(df)


items = apriori(df, min_support=0.5, use_colnames=True)

print(items)
test_df = df[0:4]
print(test_df)

print(support(items, test_df))

