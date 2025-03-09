# Get the list of all tissue names
gsutil ls gs://sqtl_gtex/files/testing/mixture_testing/real_data_outputs/ | \
grep '_alt_ss_EF_5_FPR_5_' | sed -E 's#.*/([^/]+)_alt_ss_EF_5_FPR_5_.*_alphas_betas.csv#\1#' | sort | uniq | \
while read tissue; do
    echo "Processing tissue: $tissue"
    gsutil cat gs://sqtl_gtex/files/testing/mixture_testing/real_data_outputs/"${tissue}"_alt_ss_EF_5_FPR_5_*_alphas_betas.csv | \
    gsutil cp - gs://sqtl_gtex/files/testing/mixture_testing/real_data_outputs/"${tissue}"_alt_ss_EF_5_FPR_5_all_alphas_betas.csv
done

