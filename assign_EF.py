#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import scipy
from scipy.stats import beta
import argparse


# In[2]:


parser = argparse.ArgumentParser()
parser.add_argument("--tissue_file", type=str)
parser.add_argument("-o", "--output_file", help="Directs the output to a name of your choice")


args = parser.parse_args()

tissue_file=args.tissue_file


# In[5]:


def calc_EF(EF_percent, a, b):
    
    
    x = np.linspace(0,1, 1000)
    
    y =  1-beta.cdf(x, a, b)

    
    EF = y[x>=EF_percent].max()
    return EF


# In[6]:



optimizer_results = pd.read_csv(tissue_file, low_memory=False)

# In[10]:


EF_percent=0.01
new_EF_1= optimizer_results.apply(lambda x: calc_EF(EF_percent, x['alpha_prior'], x['beta_prior']), axis=1)
optimizer_results=optimizer_results.assign(EF_1=new_EF_1)

EF_percent=0.05
new_EF_5= optimizer_results.apply(lambda x: calc_EF(EF_percent, x['alpha_prior'], x['beta_prior']), axis=1)
optimizer_results=optimizer_results.assign(EF_5=new_EF_5)

EF_percent=0.1
new_EF_10= optimizer_results.apply(lambda x: calc_EF(EF_percent, x['alpha_prior'], x['beta_prior']), axis=1)
optimizer_results=optimizer_results.assign(EF_10=new_EF_10)

EF_percent=0.2
new_EF_20= optimizer_results.apply(lambda x: calc_EF(EF_percent, x['alpha_prior'], x['beta_prior']), axis=1)
optimizer_results=optimizer_results.assign(EF_20=new_EF_20)


#convergence, filter for background and nonbackground set

converged=optimizer_results.function_output_emp_bayes.apply(lambda x: x.split('\n  ')[1]=='success: True')


optimizer_results_converged=optimizer_results[converged]


passed_EF_5_threshold=optimizer_results_converged.EF_5 >= 1/optimizer_results.number_people_in_sample.iloc[1]



optimizer_results_converged_converged_passed_threshold=optimizer_results_converged[passed_EF_5_threshold]


# In[14]:


optimizer_results_converged=optimizer_results_converged.assign(passed_min_threshold_in_tissue=passed_EF_5_threshold)


# In[15]:


#write output df
with open(args.output_file, 'w') as output_file:
    optimizer_results_converged.to_csv(output_file)
    


# In[ ]:




