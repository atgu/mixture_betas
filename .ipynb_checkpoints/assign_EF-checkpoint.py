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





def mixture_calc_EF(EF_percent, params):
  
  
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
        psi = np.linspace(0, 1, 100)
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
    
    for f in [files]:
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
    # List of columns to convert
    cols_to_convert = ['params', 'LogLikelihood',
           'LogLikelihood_triple', 'LogLikelihood_double', 'params_triple',
           'params_double']
    

    all_EFs_df[cols_to_convert] = all_EFs_df[cols_to_convert].applymap(lambda s: np.fromstring(s.strip("[]"), sep=" ") if isinstance(s, str) else np.nan)


    return all_EFs_df


files = [
    input_directory + 'SE_' + tissue_name + '.csv.gz',
    input_directory +  tissue_name + '.csv.gz'
]
# In[15]:
SEs=get_EFs_dfs(files[0])
SEs=SEs

alt_ss=get_EFs_dfs(files[1])

df_w_params_uncleaned = pd.concat([SEs, alt_ss])

df_w_params=df_w_params_uncleaned[df_w_params_uncleaned.total_reads_spanning_all_junctions!='total_reads_spanning_all_junctions']


def calc_EF(df):
    EF_1 = df.params_triple.apply(lambda x: mixture_calc_EF(0.01,x))
    
    df=df.assign(EF_1=EF_1)
    
    EF_5 = df.params_triple.apply(lambda x: mixture_calc_EF(0.05,x))
    
    df=df.assign(EF_5=EF_5)
    
    EF_10 = df.params_triple.apply(lambda x: mixture_calc_EF(0.1,x))
    
    df=df.assign(EF_10=EF_10)
    
    EF_20 = df.params_triple.apply(lambda x: mixture_calc_EF(0.2,x))
    
    df=df.assign(EF_20=EF_20)

    return df

df_w_params_w_EFs = calc_EF(df_w_params)

#dict_of_params=dict(zip(df_w_params.index, df_w_params.params_triple))

#params_as_matrix = np.array(list(dict_of_params.values()))  # Shape (rows, 9)

#matrix = np.zeros((len(dict_of_params), 4))  # Ensure correct shape


#for j, percent in enumerate([0.01, 0.05, 0.1, 0.2]):
    
    #percent_array = np.tile(percent, (len(dict_of_params.keys()), 1))


    #matrix[:, j] = np.array([mixture_calc_EF(row, params) for row, params in zip(percent_array, params_as_matrix)])

#EFs=pd.DataFrame({'EF_1':matrix[:,0],
           # 'EF_5':matrix[:,1],
            #'EF_10':matrix[:,2],
            #'EF_20':matrix[:,3]})


#df_w_params_w_EFs=pd.concat([df_w_params, EFs], axis=1)




#optimizer_results_converged=optimizer_results_converged.assign(passed_min_threshold_in_tissue=passed_EF_5_threshold)



#SEs=SEs[SEs.params_triple.apply(lambda x: len(str(x))>50 ) ]
#print('doesnt work params')
#print(sum(SEs.params_triple.apply(lambda x: len(x))!=9 ))
#alt_ss=alt_ss[alt_ss.params_triple.apply(lambda x: len(str(x))>50 ) ]

#print(alt_ss.params_triple)
     
#SEs_EFs = assign_EF_matrcies(SEs)

#alt_ss_EFs = assign_EF_matrcies(alt_ss)
#write output df
#pd.concat([SEs_EFs,alt_ss_EFs]).
df_w_params_w_EFs.to_csv(args.output_file, compression='gzip', index=False)

    


# In[ ]:




