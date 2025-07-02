# mixture_betas
This repo contains code to optimize parameters for single and mixture of betas in run alt ss w mixed effects.py

## Step-by-Step Instructions
The input files are available on GTEx website. 

https://gtexportal.org/home/downloads/adult-gtex/qtl > https://gtexportal.org/home/downloads/adult-gtex/qtl#qtl-gtex_analysis_v8-single-tissue_cis-qtl_data > GTEx_Analysis_v8_sQTL_phenotype_matrices.tar > *_perind_numers.counts.gz files. 

For ease of performing a tutorial, we have included a sample file (Brain_Cerebellar_Hemisphere) and also include the final file we used after fitting. It is the same name as on the GTEx swebsite

### Fit EFs
1. to fit parameters for potential NVEs in a tissue, you run:

```
python3 run_optimizer_alt_ss_w_mixed_effects.py --tissue "Brain_Cerebellar_Hemisphere" --filter_file alt_ss_for_EF_algorithm.csv.gz --tissue_file Brain_Cerebellar_Hemisphere_perind_numers.counts.gz --output_file output_alt_ss.csv
```
Note that you can do the same thing with skipped exons, you just need to run 

```
python3 run_optimizer_SE_w_mixed_effects.py --tissue "Brain_Cerebellar_Hemisphere" --filter_file skipped_exons_w_gene_names_and_assigned_abc_introns.csv.gz --tissue_file Brain_Cerebellar_Hemisphere_perind_numers.counts.gz --output_file output_SE.csv
```

2. then assign EFs:

```
python3 assign_EF.py --tissue_name "Brain_Cerebellar_Hemisphere" --input_directory directory_contains_step_1_files --tissue_file Brain_Cerebellar_Hemisphere_perind_numers.counts.gz --output_file Brain_Cerebellar_Hemisphere_EF_values.csv

```


## note on other files in this repo
### dockerfile
This repo also contains the dockerfile to be able to run this code, if you use docker. And the example notebook that sent the dockerimage of mixture_betas to a cluster that then runs the code
01_hail_run_pipeline.ipynb


### filter splicing events
to filter splicing events, we used the following script: post_process_EF_fits.ipynb. This is still being adapted for others to use. 
