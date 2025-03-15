#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os 
import pandas as pd
import numpy as np
import scipy
import fsspec
import argparse
import gcsfs
from scipy.stats import beta

# In[2]:


parser = argparse.ArgumentParser()
parser.add_argument("--input_directory", type=str)
parser.add_argument("--tissue_name", type=str)
parser.add_argument("-o", "--output_file", help="Directs the output to a name of your choice")


args = parser.parse_args()

tissue_name=args.tissue_name

input_directory=args.input_directory



# In[5]:



def mixture_calc_EF(EF_as_percent, params):
    """
    Compute the exceedance fraction (EF) for a mixture of Beta distributions.
    
    The function supports the following parameter orderings:
    
    1. Only distribution parameters:
         [a1, a2, ..., an,  b1, b2, ..., bn]
       (Weights are assumed equal: w_i = 1/n.)
       
    2. Fully specified weights:
         [a1, a2, ..., an,  b1, b2, ..., bn,  w1, w2, ..., wn]
       
    3. Partially specified weights (n>=2):
         [a1, a2, ..., an,  b1, b2, ..., bn,  w1, w2, ..., w_{n-1}]
       In this case, the last weight is computed as:
         w_n = 1 - (w1 + ... + w_{n-1}).
    
    The function computes the mixture cumulative distribution function (CDF) 
    on a fixed grid psi in [0,1] with 1000 points, then calculates the density 
    as:
    
         density = 1 - (weighted sum of Beta CDFs)
         
    Finally, it returns the maximum density for psi values greater than or equal 
    to EF_percent.
    
    Parameters:
        EF_percent (float): The psi threshold (between 0 and 1) at which to compute EF.
        params (list or array): The parameter vector following one of the conventions above.
    
    Returns:
        EF (float): The computed exceedance fraction.
    """
    EF_percent=EF_as_percent/100
    L = len(params)
    
    # --- Determine the parameter ordering and extract a's, b's, and weights ---
    # Case 1: Only distribution parameters (2*n values)
    if L % 2 == 0 and (L % 3 != 0 and (L + 1) % 3 != 0):
        n_components = L // 2
        a_params = params[0:n_components]
        b_params = params[n_components:2*n_components]
        w_params = [1.0 / n_components] * n_components  # equal weights
    # Case 2: Fully specified weights (3*n values)
    elif L % 3 == 0:
        n_components = L // 3
        a_params = params[0:n_components]
        b_params = params[n_components:2*n_components]
        w_params = params[2*n_components:3*n_components]
    # Case 3: Partially specified weights (3*n - 1 values)
    elif (L + 1) % 3 == 0:
        n_components = (L + 1) // 3
        a_params = params[0:n_components]
        b_params = params[n_components:2*n_components]
        w_params = list(params[2*n_components:])  # first n-1 weights
        w_last = 1 - sum(w_params)
        w_params.append(w_last)
    else:
        raise ValueError("Parameter list length not recognized. It should be either 2*n, 3*n, or 3*n-1 in length.")
    
    # --- Evaluate the mixture CDF on a fixed grid ---
    psi = np.linspace(0, 1, 1000)
    mixture_cdf = np.zeros_like(psi, dtype=float)
    
    for i in range(n_components):
        mixture_cdf += w_params[i] * beta.cdf(psi, a_params[i], b_params[i])
    
    # --- Compute density as 1 - (mixture CDF) ---
    density = 1 - mixture_cdf
    
    # --- Find the maximum density for psi >= EF_percent ---
    mask = psi >= EF_percent
    if not np.any(mask):
        raise ValueError("No psi values are greater than or equal to EF_percent.")
    EF = density[mask].max()
    
    return EF
    
    
    



 


def import_params_file(files):
    all_EFs_df = pd.DataFrame()
    
    #files = os.listdir(output_dir)
    
    for f in files:
        try:
            #if 'EF_' in f and 'alt' in f.split('EF_')[1].split('_')[0]:
            
            df_chunks = pd.read_csv(f, low_memory=False, chunksize=100000, compression='gzip')
            df = pd.concat(list(df_chunks), ignore_index=True)  # Convert chunks to a single DataFrame
                
                # Extract attributes from the filename
            tissue = f.split('mixture_')[1].split('.csv.gz')[0] if 'mixture_' in f else 'unknown'
            splicing_event = f.split('EF_')[1].split('_')[0] if 'EF_' in f else 'unknown'

            df = df[df.pvalue_dip != 'pvalue_dip']  # Remove potential header duplication
            df = df.assign(tissue=tissue, splicing_event=splicing_event)
            df['pvalue_dip'] = pd.to_numeric(df['pvalue_dip'], errors='coerce')  
            df['pvalue_LR'] = pd.to_numeric(df['pvalue_LR'], errors='coerce')  
                # Convert safely
                
            all_EFs_df = pd.concat([all_EFs_df, df], ignore_index=True)
        
        except Exception as e:
            print(f"Error processing {f}: {e}")  # Debugging log
    
    return all_EFs_df

def get_EFs_dfs(file):
    all_EFs_df = import_params_file(file)
    print(all_EFs_df.columns)
    # List of columns to convert
    cols_to_convert = ['params', 'LogLikelihood', 'LogLikelihood_triple', 'LogLikelihood_double', 'params_triple', 'params_double']
    

    all_EFs_df[cols_to_convert] = all_EFs_df[cols_to_convert].applymap(lambda s: np.fromstring(s.strip("[]"), sep=" ") if isinstance(s, str) else np.nan)


    return all_EFs_df



files = [
    input_directory + 'all_EF_SEs_alphas_and_betas_from_mixture_' + tissue_name + '.csv.gz',
    input_directory + 'all_EF_alt_ss_alphas_and_betas_from_mixture_' + tissue_name + '.csv.gz'
]

df_w_params=get_EFs_dfs(files)

dict_of_params=dict(zip(df_w_params.index, df_w_params.params_triple))

params_as_matrix = np.array(list(dict_of_params.values()))  # Shape (rows, 9)

matrix = np.zeros((len(dict_of_params), 4))  # Ensure correct shape


for j, percent in enumerate([1, 5, 10, 20]):
    
    percent_array = np.tile(percent, (len(dict_of_params.keys()), 1))


    matrix[:, j] = np.array([mixture_calc_EF(row, params) for row, params in zip(percent_array, params_as_matrix)])

EFs=pd.DataFrame({'EF_1':matrix[:,0],
            'EF_5':matrix[:,1],
            'EF_10':matrix[:,2],
            'EF_20':matrix[:,3]})


df_w_params_w_EFs=pd.concat([df_w_params, EFs], axis=1)


#optimizer_results_converged=optimizer_results_converged.assign(passed_min_threshold_in_tissue=passed_EF_5_threshold)


# In[15]:


#write output df
df_w_params_w_EFs.to_csv(args.output_file, compression='gzip', index=False)

    


# In[ ]:




