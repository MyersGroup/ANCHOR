# ANCHOR 
An approach leveraging segments of distinct ancestries within individuals to estimate similarity in underlying causal effect sizes between two groups

We recommend user to use "HAPMIX" to infer the local ancestry. If user already used "HAPMIX" and generated the output, you can use the R script: "read_imp_hapmix_16prob_P4_s1.r" to genrate ancestry specific genotypes (pop1 and pop2, pop1 should be the population on which the GWAS (effect size) was conducted). User need to provide the same arguments as what user specified in the "HAPMIX" configuration file.

Markup : * ADMIXINDFIE
		 * OUTDIR (directory of hapmix output)
		 * ADMIXPOP (hapmix output file prefix)
		 * HAPMIX_MODE (diploid as default)
		 * output.dir (store the output ancestry specific genotype and local ancestry files)

Given the above files, user can run the function "read.imp4" either within an R session by sourcing the script as follows:

```r
source('read_imp_hapmix_16prob_P4_s1.r')
read.imp4<-function(ADMIXINDFILE,OUTDIR,ADMIXPOP,HAPMIX_DATADIR,HAPMIX_MODE,output.dir,mcc=16)
``` 

Or use can run the script from the commandlines as follows:

```bash
Rscript AdmixIndFile HapmixOutDir HapmixDataDir AdmixPop Diploid test_dir 12
```

Result will be generated to the directory "test_dir"

We also highly recommend user to "mean-center" the (ancestry) genotypes condition on local ancestry before conducting any further analysis, as we noticed that the non-mean-centered two ancestry genotypes are not independent and the effect size of PGS based on non-mean-centered genotypes will be biased by the ancestry (e.g. reduced estimated values for European effect size in individuals with high African ancestry.

If User has already used the "mean-centered" score or genotype,please ignore the following section.

## (Mean-centered) Ancestry PGS construction

As illustrated in our paper, we brought up an mean-centered ancestry PGS construction methods to deconvolute ordinary PGS into two independent ancestry PGSs (e.g. European and African PGS for European-African admixed populations).

Conventional approaches to construct such Ancestry PGSs will often lead to correlation between the PGS, and therefore make the results not robust. Our method, however, can correct the potential collinearity between the ancestry PGSs by introducing so-called *Mean-centering* technique to remove the correlation between the ancestry PGSs. 

To construct the mean-centered PGS, please use the function built in the script "generate_mean_centered_PGS.r". 

We would expect the user to run "hapmix" (diploid mode) (detail can be seen at "hapmix" paper) before using this script, and our script expect output from "hapmix"

Suppose user generated "hapmix" output at "/home/userA/hapmix_test" directory. "hapmix" sample id file is "/home/userA/data/sample_ind", "HAPMIX_MODE" for "hapmix" is "deploid" and admixed population is "ADMIXPOP"="EU_AF", then user can first open an R session and source the script as follows:

```r
source('generate_mean_centered_PGS.r')
``` 

Then user can run the command as follow to generate some important intermediate fileat specified folder "/home/userA/hapmix_PGS". User can also specify the number of cores to be used to speed up the job

```r
ri<-read.imp4("/home/userA/hapmix_test/data/sample_ind","/home/userA/hapmix_test/out/","EU_AF","diploid","/home/userA/hapmix_PGS",10)
```

Afterwards, user can call function "estimate_f" to generate the effect allele frequency between two populations:

```r
fq<-estimate_f()
```

To "mean-centered" genotypes (so as to PGS), the user should run the following function:

```r
mcg<-mean.center.geno()
```

By providing a table including the effect size estimation from other GWAS analysis, in which there must be a column called "SNP" and "BETA", user can generate their own ancestry PGS as follows():

```r
my.pgs<-cal.pgs('height_beta.rds')
```

Please bear in mind that we expect the input table should be in ".rds" format. If user's original file is plain text table, user can read it into R and use "saveRDS" function to convert the format of the original table.
