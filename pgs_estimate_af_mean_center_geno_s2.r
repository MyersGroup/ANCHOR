library(parallel)

cc<-commandArgs(R)

if(length(cc)==0){
	print("We assume you would run this script within an R session")
	print("Please provide the following files:")
	print("ADMIXINDFILE=")
	print("out.geno.dir=")
	print("output.dir=")
	print("betafile=")
	print("pop1.genofile=")
}
else{
	ADMIXINDFILE<-cc[1]#'batch_1_38/data/sample_ind'
	out.geno.dir<-cc[2]#'test_s1'#this should be the output.dir in the function "read.imp4" defined in script "read_imp_hapmix_16prob_P4_s1.r"
	output.dir<-cc[3]#'test_s2'#output directory to store the output
	betafile<-cc[4]
	pop1.genofile<-cc[5]
}

pgs.anno<-function(){
	pgsanno<-readRDS(paste0(out.geno.dir,'/SNP_annot.rds'))
}
if(!exists('s.anno')){
	s.anno<-pgs.anno()
}

if(!dir.exists(output.dir)){
	dir.create(output.dir)
}

if(!exists('AFN')){
	AFN<-basename(ADMIXINDFILE)
}

estimate_f<-function(mcc=12){
	geno.info<-readRDS(paste0(out.geno.dir,'/',AFN,'.rawgeno.rds'))
	anc.A11<-readRDS(paste0(out.geno.dir,'/',AFN,'.A11.rds'))
	anc.A12<-readRDS(paste0(out.geno.dir,'/',AFN,'.A12.rds'))
	geno.info<-as.matrix(geno.info)
	anc.A11<-as.matrix(anc.A11)
	anc.A12<-as.matrix(anc.A12)
	geno.info[geno.info>2]=2
	nn<-nrow(anc.A11)
	anc.sum<-anc.A11+anc.A12
	anc.sum[anc.sum<1]=1
	anc.A11<-anc.A11/anc.sum
	anc.A12<-anc.A12/anc.sum
	anc.x<-2*anc.A11+anc.A12
	np<-mclapply(1:nn,function(i)coefficients(lm(geno.info[i,]~anc.x[i,])),mc.cores=mcc)
	np.dt<-do.call(rbind,np)
	f.pop1<-np.dt[,2]+np.dt[,1]/2
	f.pop2<-np.dt[,1]/2
	f.pop1[f.pop1<0]=0
	f.pop1[f.pop1>1]=1
	f.pop2[f.pop2<0]=0
	f.pop2[f.pop2>1]=1
	freq.dt<-s.anno
	freq.dt$pop1.freq<-f.pop1
	freq.dt$pop2.freq<-f.pop2
	freq.dt
}

if(!exists('fq')){
	if(!file.exists(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))){
		fq<-estimate_f()
		saveRDS(fq,paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))
	}
	else{
		fq<-readRDS(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))
	}
}
maf.thres=0.01

if(!exists('filter.ids')){
	pop1.maf<-pmin(fq$pop1.freq,1-fq$pop1.freq)
	pop2.maf<-pmin(fq$pop2.freq,1-fq$pop2.freq)
	filter.ids<-pop1.maf>maf.thres & pop2.maf>maf.thres
}	

if(!exists('fq.maf')){
	if(!file.exists(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq_maf.rds'))){
		fq.maf<-fq[filter.ids,]
		saveRDS(fq.maf,paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq_maf.rds'))
	}
	else{
		fq.maf<-readRDS(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq_maf.rds'))
	}
}

mean.center.geno<-function(){
	#AFN<-basename(ADMIXINDFILE)
	geno.info<-readRDS(paste0(out.geno.dir,'/',AFN,'.rawgeno.rds'))
	anc.A11<-readRDS(paste0(out.geno.dir,'/',AFN,'.A11.rds'))
	anc.A12<-readRDS(paste0(out.geno.dir,'/',AFN,'.A12.rds'))
	geno.info<-as.matrix(geno.info)
	anc.A11<-as.matrix(anc.A11)
	anc.A12<-as.matrix(anc.A12)
	geno.info[geno.info>2]=2
	nn<-nrow(anc.A11)
	anc.sum<-anc.A11+anc.A12
	anc.sum[anc.sum<1]=1
	anc.A11<-anc.A11/anc.sum
	anc.A12<-anc.A12/anc.sum
	anc.A22<-1-anc.A11-anc.A12
	ANC2<-anc.A12*0.5+anc.A22
	ANC2.global<-colMeans(ANC2)[filter.ids]#estimate global ancestry for POP2
	saveRDS(ANC2.global,paste0(output.dir,'/global_ANC_pop2.rds'))

	geno.P11<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop11.rds')))
	geno.P12<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop12.rds')))
	geno.P21<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop21.rds')))
	geno.P22<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop22.rds')))
		
	geno.P1<-geno.P11+geno.P12-(2*fq$pop1.freq*anc.A11+fq$pop1.freq*anc.A12)
	geno.P2<-geno.P22+geno.P21-(2*fq$pop2.freq*anc.A22+fq$pop2.freq*anc.A12)
	
	saveRDS(geno.P1[filter.ids,],paste0(out.geno.dir,'/',AFN,'.pop1.rds'))
	saveRDS(geno.P2[filter.ids,],paste0(out.geno.dir,'/',AFN,'.pop2.rds'))
}

align.beta<-function(beta.df){#beta.df should be a data frame
	beta.ids.1<-paste0(beta.df$chr,':',beta.df$pos_grch37,':',beta.df$other_allele,':',beta.df$effect_allele)
	beta.ids.2<-paste0(beta.df$chr,':',beta.df$pos_grch37,':',beta.df$effect_allele,':',beta.df$other_allele)
	fq.maf.snp<-paste0(fq.maf[[2]],':',fq.maf[[4]],':',fq.maf[[5]],':',fq.maf[[6]])
	
	beta.df[beta.ids.2%in%fq.maf.snp,'BETA']<- -beta.df[beta.ids.2%in%fq.maf.snp,'BETA']#flip the effect size if the ref allele in fq is the effect allele
	#swap the other/effect allele
	beta.df.tmp<-beta.df[beta.ids.2%in%fq.maf.snp,'other_allele']
	beta.df[beta.ids.2%in%fq.maf.snp,'other_allele']<-beta.df[beta.ids.2%in%fq.maf.snp,'effect_allele']
	beta.df[beta.ids.2%in%fq.maf.snp,'effect_allele']<-beta.df.tmp
	beta.df.match<-beta.df[beta.ids.1%in% fq.maf.snp | beta.ids.2%in% fq.maf.snp,]
	beta.df.ids<-paste0(beta.df.match$chr,':',beta.df.match$pos_grch37,':',beta.df.match$other_allele,':',beta.df.match$effect_allele)
	fq.beta.ids<-match(beta.df.ids,fq.maf.snp)
	beta.list<-list(beta.ndf=beta.df.match,fq.beta.ids=fq.beta.ids)
}

cal.admix.pgs<-function(betafile,outpgs.prefix=NA){#betafile must have 4 columns: 'chr','pos_grch37','other_allele(non-effect allele)','effect_allele', and one column called "BETA" including the BETA estimation from the internal/external source of gwas summary statistics
#User should also make sure the SNP column in the betafile has the following format: chr:pos_GRCh37:non_effect_allele:effect_allele
	betas<-readRDS(betafile)
	#fq.maf.snp<-paste0(fq.maf[[2]],':',fq.maf[[4]],':',fq.maf[[5]],':',fq.maf[[6]])
	#betas.in<-intersect(betas$SNP,fq.maf.snp)
	#beta<-betas[match(betas.in,betas$SNP),'BETA']
	#fq.ids<-match(betas.in,fq$SNP)
	fbb<-align.beta(betas)
	beta<-fbb$beta.ndf$BETA
	fq.ids<-fbb$fq.beta.ids
	g.p1<-readRDS(paste0(out.geno.dir,'/',AFN,'.pop1.rds'))[fq.ids,]
	g.p2<-readRDS(paste0(out.geno.dir,'/',AFN,'.pop2.rds'))[fq.ids,]
	pgs.p1<-t(g.p1)%*%matrix(beta,ncol=1)
	pgs.p2<-t(g.p2)%*%matrix(beta,ncol=1)
	pgs<-data.frame(pgs.pop1=pgs.p1,pgs.pop2=pgs.p2)
	if(is.na(outpgs.prefix)){
		saveRDS(pgs,paste0(output.dir,'/',basename(betafile),'.admixpgs.rds'))
	}
	else{
		saveRDS(pgs,paste0(output.dir,'/',outpgs.prefix,'.rds'))
	}
}

cal.external.pgs<-function(betafile,pop1.genofile,outpgs.prefix=NA){#betafile is the same as above function, the external population, should be pop1 (noramlly European population, where the betafile comes from. 
# Genotype used to constructe the External(pop1) pgs, however, doesn't need to be mean-centered.

#The genotype of the external(pop1) should match the genotype of the admixed population (used in function cal.admix.pgs)
#pop1.genofile should be MxP matrix where M is the number of markers (the same marker sets as the fq set, and should be stored in rds format 

	betas<-readRDS(betafile)
	#betas.in<-intersect(betas$SNP,fq$SNP)
	#beta<-betas[match(betas.in,betas$SNP),'BETA']
	#fq.ids<-match(betas.in,fq$SNP)
	fbb<-align.beta(betas)
	beta<-fbb$beta.ndf$BETA
	fq.ids<-fbb$fq.beta.ids
	
	g.pop1<-readRDS(pop1.genofile)

	pgs.pop1<-data.frame(pgs=t(g.pop1)%*%matrix(beta,ncol=1))
	if(is.na(outpgs.prefix)){
		saveRDS(pgs.pop1,paste0(output.dir,'/',basename(betafile),'.pop1pgs.rds'))
	}
	else{
		saveRDS(pgs.pop1,paste0(output.dir,'/',outpgs.prefix,'.rds'))
	}
}

if(length(cc)>0){
	mc<-mean.center.geno()
	admix.pgs<-cal.admix.pgs(betafile)
	pop1.pgs<-cal.external.pgs(betafile,pop1.genofile)
	q(save='no')
}
