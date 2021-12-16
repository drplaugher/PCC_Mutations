# This code runs the Kruskal-wallis test on the key mutation groups
# Written by Daniel Plaugher 2021
import pandas as pd
from pingouin import kruskal
# makes data frame for PCC expressions
df = pd.read_csv('PCC_expressionsRed.csv')



cols= df.columns.tolist() # converting columnn names to a list
cols.remove('Patient') # removing unneeded cols
cols.remove('Code')
cols.remove('Mutation')
columns=cols


#----------------- perform Kruskal-Wallis Test ----------------- 
KW_Results=[] #preassigns space for results
for x in columns: # takes column names
    z=kruskal(data=df, dv=x, between='Mutation')
    KW_Results.append(z)



#----------------- Plotting ----------------- 
import seaborn as sns

sns.catplot(x='Mutation', y='MDM2', data=df, kind='violin')
sns.catplot(x='Mutation', y='BAX', data=df, kind='violin')
sns.catplot(x='Mutation', y='PIK3CD', data=df, kind='violin')






