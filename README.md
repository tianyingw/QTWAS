<!-- README.md is generated from README.Rmd. Please edit that file -->
QTWAS
==========

We provide a quantile tool for investigating nonlinear gene-trait associations in Transcriptome-Wide Association Analysis (TWAS). More information can be found on the webpage: https://tianyingw.github.io/QTWAS/.

Download
------------

- You can download the **development** version of the main function "help_funcions.R" from [Github](https://github.com/tianyingw/QTWAS).
- You can download the pre-trained models for QTWAS from [Dropbox]()

Usage
-----
Suppose you have a file for GWAS summary statistics stored in a txt file with four columns: SNP for rsid, A1 for ref allele, A2 for alt allele, and Z for z score. Then, a basic tutorial is demonstrated from below.

``` r
# import required package
library(RSQLite)
source("help_functions.R")

# specify the name of your GWAS summary stats file name
phecode_long = "sumstat_PD" 
# indicate the column for SNP(rsid)
snp.col = 1 
# indicate the column name for rsid
snp.col.name = "SNP" 

# specify your tissue name from the list here
tissue.name = "Whole_Blood" 

# specify your phenotype here
phecode = "PD" 

# specify the name for your results file
special.end =paste0("V8_", phecode) 

# specify directories
add_savepval = "~/Dropbox/somewhere_to_save_p_values/"
add_GWAS = "~/Dropbox/somewhere_you_saved_gwas_sumstats/"
add_savemodel = "~/Dropbox/somewhere_you_saved_pretrained_QTWAS_models/"

# apply QTWAS pre-trained models
pval.mat = apply_QTWAS(tissue.name, snp.col, snp.col.name, special.end, phecode_long, add_savemodel, add_GWAS, add_savepval)

```

Output
-------
This will return a matrix with 11 columns:
- *qr*: unified p values per gene;
- *qr(0.05,0.35)*: quantile interval specific p value for the quantile region (0.05, 0.35);
- *qr(0.25,0.55)*: quantile interval specific p value for the quantile region (0.25, 0.55);
- *qr(0.45,0.75)*: quantile interval specific p value for the quantile region (0.45, 0.75);
- *qr(0.65,0.95)*: quantile interval specific p value for the quantile region (0.65, 0.95);
- *n_used*: the number of variants finally used for QTWAS test statistics;
- *n_model*: the number of variants in pre-trained QTWAS models;
- *Z_Q1*: Z score for the quantile region (0.05, 0.35);
- *Z_Q2*: Z score for the quantile region (0.25, 0.55);
- *Z_Q3*: Z score for the quantile region (0.45, 0.75);
- *Z_Q4*: Z score for the quantile region (0.65, 0.95);


License
-------

This R file is free and open source software, licensed under GPL (&gt;=2).
