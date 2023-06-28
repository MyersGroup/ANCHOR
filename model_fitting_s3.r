library(parallel)

#In this script, we will build a data table and fit the model. To estimate the confidence interval precisely, we used the bootstrap with 1000 resampling.

pop1.pheno.file<-''
pop1.covar.file<-''
pop1.pgs.file<-''

admix.pheno.file<-''
admix.covar.file<-''
admix.pgs.file<-''
admix.anc.file<-''

centralised.admix.pgs<-function(dt){
#this function is used to centralised the pgs by regressing out the 
	pop1.pgs.lm<-lm(pop1.pgs ~ .,data=t.dt[,-c(1,3,4)],na.action=na.exclude)
	pop1.pgs.rs<-resid(pop1.pgs.lm)

	pop2.pgs.lm<-lm(pop2.pgs ~ .,data=t.dt[,-c(1,2,4)],na.action=na.exclude)
	pop2.pgs.rs<-resid(pop2.pgs.lm)

	dt$pop1.pgs<-pop1.pgs.rs
	dt$pop2.pgs<-pop2.pgs.rs
	return(dt)
}	

centralised.pop1.pgs<-function(dt){
#this function is used to centralised the pgs by regressing out the 
	pop1.pgs.lm<-lm(pop1.pgs ~ .,data=t.dt[,-1],na.action=na.exclude)
	pop1.pgs.rs<-resid(pop1.pgs.lm)
	dt$pop1.pgs<-pop1.pgs.rs
	return(dt)
}	

create.pop1.df<-function(seednum=2023,mcc=12){
	#covar file should at least include age and sex, and would be great if the user can provide global ancestry (ACs or PCs) variables.
	#each file should have an identifier columns ID which contains the sample ID information
	#we ask for the rds file as input 
	pop1.pheno<-readRDS(pop1.pheno.file)
	pop1.covar<-readRDS(pop1.covar.file)
	pop1.pgs<-readRDS(pop1.pgs.file)
	
	pop1.all<-cbind(pop1.pheno,pop1.pgs[,names(pop1.pgs)!='ID',drop=F],pop1.covar[,names(pop1.covar)!='ID',drop=F])
	pop1.all<-centralised.pop1.pgs(pop1.all)
	
	bootstrap.file<-paste0(dirname(pop1.pheno.file),'/pop1_bootstrap_ids.rds')
	set.seed(seednum)
	sp.num<-nrow(pop1.all)
	if(!file.exists(bootstrap.file)){
		bootids<-sapply(1:1000,function(i)bids<-sample(1:sp.num,sp.num,replace=T))
		saveRDS(bootids,bootstrap.file)
	}
	else{
		bootids<-readRDS(bootstrap.file)
	}
	boot.dfs<-lapply(1:1000,function(i)pop1.all[bootids[,i],])
	pop1.df.list<-list(real=pop1.all,boots=boot.dfs)
}

create.admix.df<-function(seednum=2023,mcc=12){
	#covar file should at least include age and sex, and would be great if the user can provide global ancestry (ACs or PCs) variables.
	#each file should have an identifier columns ID which contains the sample ID information
	#we ask for the rds file as input 
	admix.pheno<-readRDS(admix.pheno.file)
	admix.covar<-readRDS(admix.covar.file)
	admix.pgs<-readRDS(admix.pgs.file)
	admix.anc<-readRDS(admix.anc.file)

	admix.all<-cbind(admix.pheno,admix.pgs[,names(admix.pgs)!='ID',drop=F],admix.anc[,names(admix.pgs)!='ID',drop=F],admix.covar[,names(admix.covar)!='ID',drop=F])
	admix.all<-centralised.pop1.pgs(admix.all)

	bootstrap.file<-paste0(dirname(admix.pheno.file),'/admix_bootstrap_ids.rds')
	set.seed(seednum)
	sp.num<-nrow(admix.all)
	if(!file.exists(bootstrap.file)){
		bootids<-sapply(1:1000,function(i)bids<-sample(1:sp.num,sp.num,replace=T))
		saveRDS(bootids,bootstrap.file)
	}
	else{
		bootids<-readRDS(bootstrap.file)
	}
	boot.dfs<-lapply(1:1000,function(i)admix.all[bootids[,i],])
	admix.df.list<-list(real=admix.all,boots=boot.dfs)
}

pop1.coeff<-function(seednum=2023,mcc=12){
	cpd<-create.pop1.df(seednum,mcc)
	pop1.coeff<-function(dt){
		ss<-summary(lm(pheno ~ .,data=dt))
		sq<-ss$coefficients['pop1.pgs',1:2]
		names(sq)<-c('pop1.pgs.beta','pop1.pgs.sd')
		sq
	}
	pop1.real<-pop1.coeff(cpd$real)
	pop1.boots<-do.call(rbind,mclapply(cpd$boots,pop1.coeff,mc.cores=mcc))
	pop1.coef<-list(real=pop1.real,boots=pop1.boots)
}

admix.coeff<-function(seednum=2023,mcc=12){
	cad<-create.admix.df(seednum,mcc)
	admix.coeff<-function(dt){
		ss<-summary(lm(pheno ~ .,data=dt))
		sq<-as.vector(ss$coefficients[c('pop1.pgs','pop2.pgs'),1:2])
		names(sq)<-c('pop1.pgs.beta','pop2.pgs.beta','pop2.pgs.sd','pop2.pgs.sd')
		sq
	}
	admix.real<-admix.coeff(cad$real)
	admix.boots<-do.call(rbind,mclapply(cad$boots,admix.coeff,mc.cores=mcc))
	admix.coef<-list(real=admix.real,boots=admix.boots)
}

if(!exists('beta.obs')){
	beta.obs<-pop1.coeff()
}

if(!exists('beta.pop1')){
	beta.pop1<-admix.coeff()
}

estimate.ratio<-function(mcc=12){
	obs.real<-beta.obs$real
	obs.boot<-beta.obs$boots
	pop1.real<-beta.pop1$real
	pop1.boot<-beta.pop1$boots

	rho.real<-pop1.real[['pop1.pgs.beta']]/obs.real[['pop1.pgs.beta']]

	rho.boots<-unlist(mclapply(1:length(obs.boot),function(i){
		rho.i<-pop1.boot[i,'pop1.pgs.beta']/obs.boot[i,'pop1.pgs.beta']
	},mc.cores=mcc))
	
	rho.5pctCI<-quantile(rho.boots,c(0.025,0.975))
	rho.info<-c(rho.real,rho.5pctCI)
	names(rho.info)<-c('real ratio','2.5pctCI','97.5pctCI')
	rho.info
}
