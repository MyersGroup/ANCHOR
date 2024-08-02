library(parallel)

cc<-commandArgs(T)
#In this script, we will build a data table and fit the model. To estimate the confidence interval precisely, we used the bootstrap with boots.num=1000(default) resampling.

if(length(cc)==0){
	print("We assume you would run this script within an R session")
	print("Please provide the following files:")
	print("pop1.pheno.file=")
	print("pop1.covar.file=")
	print("pop1.pgs.file=")
	print("admix.pheno.file=")
	print("admix.covar.file=")
	print("admix.pgs.file=")
	print("admix.anc.file=")
	print("boots.num=")
}else{
	pop1.pheno.file<-cc[1]
	pop1.covar.file<-cc[2]
	pop1.pgs.file<-cc[3]

	admix.pheno.file<-cc[4]
	admix.covar.file<-cc[5]
	admix.pgs.file<-cc[6]
	boots.num<-as.integer(cc[7])#optional, default is 1000
	admix.anc.file<-cc[8] #this argument can be NA/NULL
	outfile<-cc[9] #if this is NA/NULL, then the result will be redirect to the stdout (on screen) and save to an rds object at the folder where you run the script; otherwise the result will be write to the outfile
	print(paste0("outfile is:", outfile))
}	

if(!exists('boots.num') || is.na(boots.num)){
	boots.num=1000
}

centralised.admix.pgs<-function(t.dt,ifpop2.anc=F){
#this function is used to centralised the pgs by regressing out the 
	if(ifpop2.anc){
		rm.pop1.cols<-c(1,3,4)
		rm.pop2.cols<-c(1,2,4)
	}
	else{
		rm.pop1.cols<-c(1,3)
		rm.pop2.cols<-c(1,2)
	}
	pop1.pgs.lm<-lm(pop1.pgs ~ .,data=t.dt[,-rm.pop1.cols],na.action=na.exclude)
	pop1.pgs.rs<-resid(pop1.pgs.lm)

	pop2.pgs.lm<-lm(pop2.pgs ~ .,data=t.dt[,-rm.pop2.cols],na.action=na.exclude)
	pop2.pgs.rs<-resid(pop2.pgs.lm)

	t.dt$pop1.pgs<-pop1.pgs.rs
	t.dt$pop2.pgs<-pop2.pgs.rs
	return(t.dt)
}	

centralised.pop1.pgs<-function(t.dt){
#this function is used to centralised the pgs by regressing out the 
	pop1.pgs.lm<-lm(pop1.pgs ~ .,data=t.dt[,-1],na.action=na.exclude)
	pop1.pgs.rs<-resid(pop1.pgs.lm)
	t.dt$pop1.pgs<-pop1.pgs.rs
	return(t.dt)
}	

create.pop1.df<-function(seednum=2023,mcc=12){
	#covar file should at least include age and sex, and would be great if the user can provide global ancestry (ACs or PCs) variables.
	#each file should have an identifier columns ID which contains the sample ID information
	#we ask for the rds file as input 
	if(!exists('pop1.pheno.file')){
		stop("pop1.pheno.file hasn't been specified, please sepecify the pop1.pheno.file")
	}
	if(!exists('pop1.covar.file')){
		stop("pop1.covar.file hasn't been specified, please sepecify the pop1.covar.file")
	}
	if(!exists('pop1.pgs.file')){
		stop("pop1.pgs.file hasn't been specified, please sepecify the pop1.pgs.file")
	}
	
	pop1.pheno<-readRDS(pop1.pheno.file)
	pop1.covar<-readRDS(pop1.covar.file)
	pop1.pgs<-readRDS(pop1.pgs.file)

	if(!identical(pop1.pheno$ID,pop1.covar$ID) && !identical(pop1.pheno$ID,pop1.pgs$ID)){
		print('Warning: the ID order between files are not identical!')
		print('Will try to align the order')
		ID.c<-intersect(pop1.pheno$ID,intersect(pop1.covar$ID,pop1.pgs$ID))
		pop1.pheno<-pop1.pheno[match(ID.c,pop1.pheno$ID),]
		pop1.covar<-pop1.covar[match(ID.c,pop1.covar$ID),]
		pop1.pgs<-pop1.pgs[match(ID.c,pop1.pgs$ID),]
	}
	
	pop1.all<-cbind(pop1.pheno,pop1.pgs[,names(pop1.pgs)!='ID',drop=F],pop1.covar[,names(pop1.covar)!='ID',drop=F])
	pop1.all$ID<-NULL
	names(pop1.all)[1:2]<-c('pheno','pop1.pgs')
	pop1.all<-centralised.pop1.pgs(pop1.all)
	
	bootstrap.file<-paste0(dirname(pop1.pheno.file),'/pop1_bootstrap_ids.rds')
	set.seed(seednum)
	sp.num<-nrow(pop1.all)
	if(!file.exists(bootstrap.file)){
		bootids<-sapply(1:boots.num,function(i)bids<-sample(1:sp.num,sp.num,replace=T))
		saveRDS(bootids,bootstrap.file)
	}
	else{
		bootids<-readRDS(bootstrap.file)
	}
	boot.dfs<-lapply(1:boots.num,function(i)pop1.all[bootids[,i],])
	pop1.df.list<-list(real=pop1.all,boots=boot.dfs)
}

create.admix.df<-function(ifpop2.anc=F,seednum=2023,mcc=12){
	#covar file should at least include age and sex, and would be great if the user can provide global ancestry (ACs or PCs) variables.
	#each file should have an identifier columns ID which contains the sample ID information
	#we ask for the rds file as input 
	if(!exists('admix.pheno.file')){
		stop("admix.pheno.file hasn't been specified, please sepecify the admix.pheno.file")
	}
	if(!exists('admix.covar.file')){
		stop("admix.covar.file hasn't been specified, please sepecify the admix.covar.file")
	}
	if(!exists('admix.pgs.file')){
		stop("admix.pgs.file hasn't been specified, please sepecify the admix.pgs.file")
	}
	admix.pheno<-readRDS(admix.pheno.file)
	admix.covar<-readRDS(admix.covar.file)
	admix.pgs<-readRDS(admix.pgs.file)
	if(!identical(admix.pheno$ID,admix.covar$ID) && !identical(admix.pheno$ID,admix.pgs$ID)){
		print('Warning: the ID order between files are not identical!')
		print('Will try to align the order')
		ID.c<-intersect(admix.pheno$ID,intersect(admix.covar$ID,admix.pgs$ID))
		admix.pheno<-admix.pheno[match(ID.c,admix.pheno$ID),]
		admix.covar<-admix.covar[match(ID.c,admix.covar$ID),]
		admix.pgs<-admix.pgs[match(ID.c,admix.pgs$ID),]
	}
	if(ifpop2.anc){
		if(!exists('admix.anc.file')){
			stop("admix.anc.file hasn't been specified, please sepecify the admix.anc.file")
		}
		admix.anc<-readRDS(admix.anc.file)
		admix.all<-cbind(admix.pheno,admix.pgs[,names(admix.pgs)!='ID',drop=F],admix.anc[,names(admix.anc)!='ID',drop=F],admix.covar[,names(admix.covar)!='ID',drop=F])
		admix.all$ID<-NULL
		names(admix.all)[1:4]<-c('pheno','pop1.pgs','pop2.pgs','pop2.ancestry')
		admix.all<-centralised.admix.pgs(admix.all,ifpop2.anc=T)
	}
	else{
		admix.all<-cbind(admix.pheno,admix.pgs[,names(admix.pgs)!='ID',drop=F],admix.covar[,names(admix.covar)!='ID',drop=F])
		admix.all$ID<-NULL
		names(admix.all)[1:3]<-c('pheno','pop1.pgs','pop2.pgs')
		admix.all<-centralised.admix.pgs(admix.all,ifpop2.anc=F)
	}
	
	bootstrap.file<-paste0(dirname(admix.pheno.file),'/admix_bootstrap_ids.rds')
	set.seed(seednum)
	sp.num<-nrow(admix.all)
	if(!file.exists(bootstrap.file)){
		bootids<-sapply(1:boots.num,function(i)bids<-sample(1:sp.num,sp.num,replace=T))
		saveRDS(bootids,bootstrap.file)
	}
	else{
		bootids<-readRDS(bootstrap.file)
	}
	boot.dfs<-lapply(1:boots.num,function(i)admix.all[bootids[,i],])
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

admix.coeff<-function(ifpop2.anc=F,seednum=2023,mcc=12){
	cad<-create.admix.df(ifpop2.anc,seednum,mcc)
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

estimate.ratio<-function(mcc=12){
	if(!exists('beta.pop1')){
		beta.pop1<-pop1.coeff()
	}

	if(!exists('beta.admix')){
		beta.admix<-admix.coeff()
	}

	obs.real<-beta.pop1$real
	obs.boot<-beta.pop1$boots
	pop1.real<-beta.admix$real
	pop1.boot<-beta.admix$boots

	rho.real<-pop1.real[['pop1.pgs.beta']]/obs.real[['pop1.pgs.beta']]

	rho.boots<-unlist(mclapply(1:nrow(obs.boot),function(i){
		rho.i<-pop1.boot[i,'pop1.pgs.beta']/obs.boot[i,'pop1.pgs.beta']
	},mc.cores=mcc))
	
	rho.5pctCI<-quantile(rho.boots,c(0.025,0.975))
	rho.info<-c(rho.real,rho.5pctCI)
	names(rho.info)<-c('real_ratio','2.5pctCI','97.5pctCI')
	rho.info
}

get.trait.matrix.central<-function(seednum=2023,mcc=12,bootstrap=NA){
	if(!exists('tgmc')){
		tgmc<<-create.admix.df(ifpop2.anc=T,seednum=seednum,mcc=mcc)
	}
	if(is.na(bootstrap)){
		tgmc.ret<-tgmc$real
	}
	else{
		tgmc.ret<-tgmc$boots[[bootstrap]]
	}
	tgmc.ret
}

fitmodel<-function(phen,R,a=seq(-1,1,0.02)){
#this function is used to fit a non-linear model by taking into account both local and global ancestry
	mat=as.matrix(R)
	mat[,1]=rowMeans(R[,1:2])
	mat[,2]=0.5*(R[,1]-R[,2])
	coeffmat=matrix(nrow=0,ncol=ncol(R))
	ss=vector(length=0)
	oldmat=mat
	for(i in 1:length(a)){
		pred=oldmat[,1]+a[i]*oldmat[,2]
		###constant term, then predictor, then interaction with ancestry, then other regressors
		mat=cbind(1,pred,pred*oldmat[,3],oldmat[,4:ncol(mat)])
		####if no interaction with ancestry needed, then fit the simplemodel
		####if only interaction with ancestry then set a=0 and fit again
		#####the simplest model is a=0 and nothing else
		#####so get 4 different fits
		aa=solve (t(mat) %*% mat)
		bb=t(mat)%*%matrix(phen,ncol=1)
		predcoeff=aa%*%bb
		coeffmat=rbind(coeffmat,as.vector(predcoeff))
		predphen=mat %*%predcoeff
		vec=phen-predphen
		ss=c(ss,sum(vec^2))
	}
	newcoeffmat=cbind(coeffmat[,1],0.5*coeffmat[,2]+0.5*a*coeffmat[,2],0.5*coeffmat[,2]-0.5*a*coeffmat[,2],coeffmat[,3]/coeffmat[,2],coeffmat[,4:ncol(coeffmat)])  #### recombine the coef for eu.score and af.score
	####if no interaction with ancestry needed, then fit the simple model
	simple=lm(phen~R[,-3])
	newcoeff=c(simple$coeff[1:3],0,simple$coeff[4:length(simple$coeff)])
	rssnointeractionwithancestry=sum(simple$resid^2)
	coeffmat=rbind(newcoeff,newcoeffmat)
	####if only interaction with ancestry then set a=0 and fit again
	pred=oldmat[,1]
	###constant term, then predictor, then interaction with ancestry, then other regressors
	mat=cbind(1,pred,pred*oldmat[,3],oldmat[,4:ncol(mat)])
	aa=solve (t(mat) %*% mat)
	bb=t(mat)%*%matrix(phen,ncol=1)
	predcoeff=aa%*%bb
	predphen=mat %*%predcoeff
	vec=phen-predphen
	rssnolocaldiffonlyancestry=sum(vec^2)
	newcoeff=c(predcoeff[1],0.5*predcoeff[2],0.5*predcoeff[2],predcoeff[3]/predcoeff[2],predcoeff[4:length(predcoeff)])###
	coeffmat=rbind(coeffmat[1,],newcoeff,coeffmat[2:nrow(coeffmat),])
	#####the simplest model is a=0 and nothing else
	simplest=lm(phen~oldmat[,-c(2,3)])
	rssnolocaldiffnointeractionwithancestry=sum(simplest$resid^2)
	newcoeff=c(simplest$coeff[1],0.5*simplest$coeff[2],0.5*simplest$coeff[2],0,simplest$coeff[3:length(simplest$coeff)])
	coeffmat=rbind(newcoeff,coeffmat)
	#####so get 4 different fits
	dfs=c(simplest$df.residual,simplest$df.residual-1,simple$df.residual,simple$df.residual-1)
	return(list(dfs=dfs,rssbasic=rssnolocaldiffnointeractionwithancestry,rssonlylocal=rssnointeractionwithancestry,rssonlyglobal=rssnolocaldiffonlyancestry,rssfullmodel=ss,alphas=a,coeffmat=coeffmat))
}

get.the.central.pv<-function(bootstrap=NA){
	gtm<-get.trait.matrix.central(bootstrap)
	gtm<-gtm[complete.cases(gtm),]
	fm<-fitmodel(gtm[,1],as.matrix(gtm[,-1]))
	rss.bas<-fm$rssbasic
	rss.locOnly<-fm$rssonlylocal
	rss.gloOnly<-fm$rssonlyglobal
	rss.fullmodel<-min(fm$rssfullmodel)
	test.vec<-c(rss.bas,rss.locOnly,rss.gloOnly,rss.fullmodel)
	RSSs<-c(rss.bas,rss.locOnly,rss.gloOnly,rss.fullmodel,var(gtm[,1],na.rm=T)*(nrow(gtm)-1))
	chi.square.stats<-c(rss.bas/fm$df[1],-(rss.locOnly-rss.bas)/rss.bas*fm$df[1],-(rss.gloOnly-rss.bas)/rss.bas*fm$df[1],(rss.bas-rss.fullmodel)/rss.bas*fm$df[1],(rss.gloOnly-rss.fullmodel)/rss.gloOnly*fm$df[3],(rss.locOnly-rss.fullmodel)/rss.locOnly*fm$df[2])
	dfs=c(0,1,1,2,1,1)#2:6
	chisq.pv<-1-pchisq(chi.square.stats,dfs)
	names(chisq.pv)<-c('Base','LocalVsBase','GlobalVsBase','BaseVsFull','GlobalVsFull','LocalVsFull')
	ak<-fm$alphas[fm$rssfullmodel==min(fm$rssfullmodel)]
	coeff.dt<-rbind(fm$coeffmat[1:3,2:4],fm$coeffmat[3+which(fm$rssfullmodel==min(fm$rssfullmodel)),2:4])
	#c0<-c(coeff.dt[4,1:2]*(1+coeff.dt[4,3]*0.0),0)
	c100<-c(coeff.dt[4,1:2]*(1+coeff.dt[4,3]*1.0),0)
	coeff.dt<-rbind(coeff.dt,c100)

	ret.coeff.mat<-rbind(fm$coeffmat[1:3,],fm$coeffmat[3+which(fm$rssfullmodel==min(fm$rssfullmodel)),])
	zm<-rep(1,ncol(ret.coeff.mat))
	zm[2:3]<-(1+coeff.dt[4,3])
	ret.coeff.mat<-rbind(ret.coeff.mat,ret.coeff.mat[4,]*zm)
	colnames(coeff.dt)<-c('eu.pgs','af.pgs','ratio')
	rownames(coeff.dt)<-c('simpest','local','global','full.100%pop1','full.100%pop2')
	attr(chi.square.stats,'df')<-dfs
	mat<-matrix(0,nrow=nrow(gtm),ncol=ncol(gtm))
	mat[,2]=rowMeans(gtm[,2:3])
	mat[,3]=0.5*(gtm[,2]-gtm[,3])
	pred=mat[,2]+ak*mat[,3]
	pred2=pred*gtm[,4]
	mat[,2]=pred
	mat[,3]=pred2
	mat[,c(1,4:ncol(gtm))]<-gtm[,c(1,4:ncol(gtm))]
	return(list(chisq.stat=chi.square.stats[2:6],coeff=coeff.dt,chisq.pv=chisq.pv[-1]))
}

central.bootstrap<-function(mcc=12){
	mk<-mclapply(1:boots.num,function(i)get.the.central.pv(bootstrap=i),mc.cores=mcc)
	mk.coeff.dt<-array(0,dim=c(5,2,boots.num))
	mk.chi.sq<-matrix(0,ncol=5,nrow=boots.num)
	colnames(mk.chi.sq)<-c('LocalVsBase','GlobalVsBase','BaseVsFull','GlobalVsFull','LocalVsFull')
	for(i in 1:boots.num){
		mk.coeff.dt[,,i]<-mk[[i]][[2]][,1:2]
		mk.chi.sq[i,]<-mk[[i]][[1]]
	}
	return(list(bootstrap.list=mk,bootstrap.chisq=mk.chi.sq,bootstrap.coeff=mk.coeff.dt))
}

generate.dt<-function(){
	real.ls<-get.the.central.pv()
	bstrp.ls<-central.bootstrap()
	real.coeff<-real.ls$coeff[1:4,]
	real.coeff<-rbind(real.coeff,real.coeff[4,]*(1+real.coeff[4,3]))
	mk.coeff.dt<-array(0,dim=c(5,7,length(bstrp.ls[[1]])))
	for(i in 1:length(bstrp.ls[[1]])){
		zu<-bstrp.ls[[1]][[i]][[2]]
		zu<-zu[1:4,]
		zu<-rbind(zu,zu[4,]*(1+zu[4,3]))
		mk.coeff.dt[,1:3,i]<-zu
		mk.coeff.dt[,4,i]<-1/(zu[4,3]+1)
		mk.coeff.dt[,5,i]<-zu[2,2]/zu[2,1]
		mk.coeff.dt[,6,i]<-zu[4,2]/zu[4,1]
		mk.coeff.dt[,7,i]<-zu[5,2]/zu[5,1]
	}
		
	real.dt<-cbind(real.coeff[,1:3],1/(real.coeff[4,3]+1),real.coeff[2,2]/real.coeff[2,1],real.coeff[4,2]/real.coeff[4,1],real.coeff[5,2]/real.coeff[5,1])
	colnames(real.dt)[4:7]<-c('ratio2','ave.ratio','pop2.ratio','pop1.ratio')
	cfv.025<-apply(mk.coeff.dt,c(1,2),function(x)quantile(x,0.025))
	cfv.975<-apply(mk.coeff.dt,c(1,2),function(x)quantile(x,0.975))
	rownames(cfv.025)<-rownames(cfv.975)<-rownames(real.dt)
	colnames(cfv.025)<-colnames(cfv.975)<-colnames(real.dt)
	plot.dt<-list(real=real.dt,ci.025=cfv.025,ci.975=cfv.975)
}

if(length(cc)>0){
	es<-estimate.ratio()
	if(is.null(outfile) || is.na(outfile)){
		print(es)
		saveRDS(es,'ratio_estimate.rds')
	}
	else{
		print(es)
		cat(paste0(paste(names(es),collapse=' '),'\n'),file=outfile)
		cat(es,file=outfile,append=T)
	}
	q(save='no')
}
