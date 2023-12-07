# This code plots Kaplan Meier curves on the key mutation groups
# Written by Daniel Plaugher 2021


#pip install lifelines
import pandas as pd
import matplotlib.pyplot as plt
import lifelines 
from lifelines import KaplanMeierFitter


# Run for combinations: TP53, Kras, T/K, S/K, T/C/K, T/S/K
df = pd.read_csv("paad_mutation_status_2018july26.csv") # read data

keep = ['Patient Number', 'Mutation', 'vital_status', 'group', 'new_days_to_death']
df=df[keep] # removes unneeded columns

# filters data fram by mutation category 
df_KRAS=df.loc[df['Mutation'] == 'KRAS']
df_TK=df.loc[df['Mutation'] == 'T/K']
df_TSK=df.loc[df['Mutation'] == 'T/S/K']
df_TCK=df.loc[df['Mutation'] == 'T/C/K']


# Overall
durations = df['new_days_to_death'] # survival data
event_observed = df['group'] #group is survival category alive (0), dead (1)


# KRAS
durations3 = df_KRAS['new_days_to_death']
event_observed3 = df_KRAS['group']

# TK
durations4 = df_TK['new_days_to_death']
event_observed4 = df_TK['group']


# TSK
durations6 = df_TSK['new_days_to_death']
event_observed6 = df_TSK['group']


# Overall
durations7 = df_TCK['new_days_to_death']
event_observed7 = df_TCK['group']



## create a kmf object
kmf = KaplanMeierFitter() 
kmf2 = KaplanMeierFitter()
kmf3 = KaplanMeierFitter()
kmf4 = KaplanMeierFitter()
kmf5 = KaplanMeierFitter()


## Fit the data into the model
kmf.fit(durations, event_observed,label='Overall')
kmf2.fit(durations3, event_observed3,label='KRAS')
kmf3.fit(durations4, event_observed4,label='T/K')
kmf4.fit(durations6, event_observed6,label='T/S/K')
kmf5.fit(durations7, event_observed7,label='T/C/K')


## Create an estimate
kmf.plot(ci_show=False) ## ci_show is meant for Confidence interval, since our data set is too tiny, thus i am not showing it.
kmf2.plot(ci_show=False) 
kmf3.plot(ci_show=False) 
kmf4.plot(ci_show=False) 
kmf5.plot(ci_show=False) 
plt.ylabel('Survival Probability')
plt.xlabel('Days')
plt.title('Kaplan Meier Plots for PDAC')
