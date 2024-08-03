# ANCHOR 
An approach leveraging segments of distinct ancestries within individuals to estimate similarity in underlying causal effect sizes between two groups

## Read output file from "HAPMIX"

We recommend user to use "HAPMIX" to infer the local ancestry. If user already used "HAPMIX" and generated the output, you can use the R script: "src/read_imp_hapmix_16prob_P4_s1.r" to generate ancestry specific genotypes (pop1 and pop2, pop1 should be the population on which the GWAS (effect size) was conducted). User need to provide the same arguments as what user specified in the "HAPMIX" configuration file.

* ADMIXINDFIE
* OUTDIR (directory of hapmix output)
* HAPMIX_DATADIR (directory of hapmix data directory)
* ADMIXPOP (hapmix output file prefix)
* HAPMIX_MODE (diploid as default)
* output.dir (store the output ancestry specific genotype and local ancestry files)

Suppose user's "hapmix" project at "/home/userA/hapmix_test" directory. "hapmix" sample id file is "/home/userA/hapmix_test/data/sample_ind";raw genotype/SNP files are at "/home/userA/hapmix_test/data/"; prefix of hapmix output file is  "AA" and HAPMIX_MODE is diploid (so the output file will be something like "AA.DIPLOID.19.20" where 19 is the chromosome number and 20 is the index of the 21st individuals (hapmix index start from 0);hapmix output path is "/home/userA/hapmix_test/output"; therefore, the variables within the script is:

```r
ADMIXINDFILE="/home/userA/hapmix_test/data/sample_ind"
OUTDIR="/home/userA/hapmix_test/output"
HAPMIX_DATADIR="/home/userA/hapmix_test/data/"
ADMIXPOP="AA"
HAPMIX_MODE="DIPLOID"
output.dir="/home/userA/hapmix_PGS"
mcc=12
```

Then user can run the command as follow to generate some important intermediate file at specified folder "/home/userA/hapmix_PGS". User can also specify the number of cores via "mcc" to speed up the job

Given the above files, user can run the function "read.imp4" either within an R session by sourcing the script as follows:

```r
source('src/read_imp_hapmix_16prob_P4_s1.r')
read.imp4<-function(ADMIXINDFILE,OUTDIR,ADMIXPOP,HAPMIX_DATADIR,HAPMIX_MODE,output.dir,mcc=16)
``` 

Or use can run the script from the command lines as follows:

```bash
Rscript src/read_imp_hapmix_16prob_P4_s1.r AdmixIndFile HapmixOutDir HapmixDataDir AdmixPop Diploid /home/userA/hapmix_PGS 12
```

Result will be generated to the directory "/home/userA/hapmix_PGS"

We also highly recommend user to "mean-center" the (ancestry) genotypes condition on local ancestry before conducting any further analysis, as we noticed that the non-mean-centered two ancestry genotypes are not independent and the effect size of PGS based on non-mean-centered genotypes will be biased by the ancestry (e.g. reduced estimated values for European effect size in individuals with high African ancestry.

If user has already used the "mean-centered" score or genotype, please ignore the following section.

## (Mean-centered) Ancestry PGS construction

As illustrated in our paper Hu et al. 2023 , we brought up an mean-centered ancestry PGS construction methods to deconvolute ordinary PGS into two independent ancestry PGSs (e.g. European and African PGS for European-African admixed populations).

Conventional approaches to construct such Ancestry PGSs will often lead to correlation between the PGS, and therefore make the results not robust. Our method, however, can correct the potential collinearity between the ancestry PGSs by introducing so-called *Mean-centering* technique to remove the correlation between the ancestry PGSs. 

In addition, we also highly recommend the user to mask the short or uncertain regions. The reasons are:

1. Our method is based on the assumption that "ancestry LD" is longer than the normal LD
2. Although the masking/unmasking results are very consitent in our simulation, we found masking help reduce the downward bias for a few of traits in the UKB real phenotypes, especially when you have some samples have almost homogeneous pop2 ancestry (e.g. >97.5% African ancestry). Centain samples should be treated as "non-admixed" samples, but sometimes "HAPMIX" may still infer some "pop1" (e.g. European) segements from those samples, which essentially will bias the estimation
3. User should avoid including such samples with extremely high proportion of "pop2" ancestry (e.g. African ancestry >95% or 97.5%), and apply both the mean-centering and masking together
4. The default masking length is 5MB (2.5MB from the central variant towards each direction), simulatin results shows this setting works very well. User can also give other setting, but please make sure the masking length is larger than 1MB (as we assume ancestry LD larger than normal LD).

To construct the mean-centered PGS, please use the function built in the script "pgs_estimate_af_mean_center_geno_s2.r". 

```r
source('src/pgs_estimate_af_mean_center_geno_s2.r')
```

User also need to provide some arguments to run the scripts, including "ADMIXINDFILE" which is the same as mentioned above, and "out.geno.dir" is the "output.dir" used in the above script (example directory:/home/userA/hapmix_PGS). User also need to provide the output directory "output.dir" and external gwas summary statistics file "betafile" to run this script  

Afterwards, a built-in function "estimate_f" will be called to generate the effect allele frequency between two populations, and the estimated allele frequency will be stored in R variable "fq":

To "mean-centered" genotypes (so as to PGS), the user should run the following function:

```r
mcg<-mean.center.geno()
```

The "betafile" is a table including the effect size estimated from other GWAS analysis, which should includes the following columns: 'chr','pos_grch37','other_allele(non-effect allele)','effect_allele', and one column called "BETA" including the BETA estimation from the internal/external source of gwas summary statistics. 

The script will try to align the gwas summary statistics and user's genotype file when user generates their own ancestry PGS as follows():

```r
admix.pgs<-cal.admix.pgs('height_beta.rds')
```

To generate the observed "PGS" (generated in the sample population as the population from which the GWAS was conducted), user need to provide 
"pop1.genofile" which is the same format as the ancestry genotype file.

Then user can generate the PGS for POP1 which has the same ancestry as the population where GWAS was conducted

```r
pop1.pgs<-cal.external.pgs(betafile,pop1.genofile)
```

Please bear in mind that we expect the input table should be in ".rds" format. If user's original file is plain text table, user can read it into R and use "saveRDS" function to convert the format of the original table.

As the script "src/read_imp_hapmix_16prob_P4_s1.r", user can also run the script "src/pgs_estimate_af_mean_center_geno_s2.r" in the command line:

```bash
Rscript src/pgs_estimate_af_mean_center_geno_s2.r AdmixIndFile /home/userA/hapmix_PGS /home/userA/mean_centered_geno "test_beta.rds" "pop1_geno.rds"
```

## Run ANCHOR to estimate effect size correlation between populations within the admixed samples

ANCHOR can leverage local ancestry information to estimate genetic correlation between samples without GWAS summary statistics from under-represent populations (POP2) in GWAS (e.g. African)."src/model_fitting_s3.r" is the main script for user to run ANCHOR method.

To run the ANCHOR method, user need to prepare the ancestry specific PGS either self-generated or using our tools. User also need to provide phenotype, any covariates (age, sex and global ancestry etc.), as well as PGS for external/POP1 where GWAS summary stats comes from, and the phenotype and covariates files (should be same to make sure the result robust). Therefore, user should provide 6 mandatory files: "pop1.pheno.file","pop1.covar.file","pop1.pgs.file","admix.pheno.file","admix.covar.file" and "admix.pgs.file". In addition, user may provide global ancestry for POP2 (if user want to fit model to predict the effect size of PGS in POP1 and POP2 with different global ancestry of POP2) file "admix.anc.file". As the confidence interval is estimated by bootstrap method, user can specify the number of bootstrap iteration, while the default value is 1000. User can also specify the output file for ratio(correlation) estimate otherwise there will be a file "ratio_estimate.rds" automatically generated at the folder where the user run this script.

Given the mandatory files provided, User can run the function "estimate.ratio" in an R session:

```r
source('src/model_fitting_s3.r')
es<-estimate.ratio()
```

or run the script from commandline:

```bash
Rscript src/model_fitting_s3.r pop1.pheno.file pop1.covar.file pop1.pgs.file admix.pheno.file admix.covar.file admix.pgs.file" 
```


## Bug report

If user comes across any problem, please email leunghom@gmail.com and we will help to resolve the problem
