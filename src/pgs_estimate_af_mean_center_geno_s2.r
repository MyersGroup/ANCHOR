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
	print("mcc=")
	print("masking.dist=")
}
else{
	ADMIXINDFILE<-cc[1]#'batch_1_38/data/sample_ind'
	out.geno.dir<-cc[2]#'test_s1'#this should be the output.dir in the function "read.imp4" defined in script "read_imp_hapmix_16prob_P4_s1.r"
	output.dir<-cc[3]#'test_s2'#output directory to store the output
	betafile<-cc[4]
	pop1.genofile<-cc[5]
	mcc<-as.integer(cc[6]])
	if(length(cc[7])==0){
		masking.dis<-2.5e6
	}
	else{
		masking.dis<-as.numeric(cc[7])
	}
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
geno.bychr.dir<-paste0(out.geno.dir,'/imp_geno/')
anc.bychr.dir<-paste0(out.geno.dir,'/imp_anc/')

estimate_f<-function(chr,mcc=12){
	geno.info<-readRDS(paste0(geno.bychr.dir,'/',AFN,'.rawgeno.chr',chr,'.rds'))
	anc.A11<-readRDS(paste0(anc.bychr.dir,'/',AFN,'.A11.chr',chr,'.rds'))
	anc.A12<-readRDS(paste0(anc.bychr.dir,'/',AFN,'.A12.chr',chr,'.rds'))
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
	freq.dt<-s.anno[s.anno[[2]]==chr,]
	freq.dt$pop1.freq<-f.pop1
	freq.dt$pop2.freq<-f.pop2
	freq.dt
}

if(!exists('fq')){
	if(!file.exists(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))){
		#fq<-estimate_f()
		fq<-do.call(rbind,lapply(1:22,estimate_f,mcc=mcc))
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
#######################
#We recommand the user to apply the masking, i.e. removing the short/uncertain ancestry segments using the function below. Based on our experiments, we thought overall 5MB should be a good threshold.
if(!exists("imp.chrs")){
	imp.chrs<-lapply(1:22,function(i)fq.maf[fq.maf[[2]]==i,])
}

get.segment.wd<-function(chr,dis=2.5e6){#get the indices of 5MB interval for each variant
	g.chr<-imp.chrs[[chr]]
	nn<-nrow(g.chr)
	cpc<-function(i){
		il=ir=NA
		if(i==1){
			il=1
		}
		if(i==nn){
			ir=nn
		}
		if(is.na(il)){#find the stretch to the left from the current position
			uz.ldf<-g.chr[i,4]-g.chr[1:(i-1),4]
			il.tmp<-which(uz.ldf < dis)
			if(length(il.tmp)==0){
				il=1
			}
			else{
				il=il.tmp[1]
			}
		}
		if(is.na(ir)){#find the stretch to the right from the current position
			uz.rdf<-g.chr[(i+1):nn,4]-g.chr[i,4]
			ir.tmp<-which(uz.rdf > dis)
			if(length(ir.tmp)==0){
				ir=nn
			}
			else{
				ir=i+ir.tmp[1]-1
			}
		}
		zp=c(il,ir)
	}
	lz<-mclapply(1:nrow(g.chr),cpc,mc.cores=mcc)
	lz.dt<-do.call(rbind,lz)
}

if(!exists('imp.anc.intv')){
	imp.anc.intv<-lapply(1:22,get.segment.wd,dis=masking.dis)
}

outpdir=paste0(out.geno.dir,'/imp_mask_region/')
if(!dir.exists(outpdir)){
	dir.create(outpdir)
}

mask.all<-function(dis.thres=masking.dis,mcc=12){
	each_chr<-function(chr){
		anc.ee<-readRDS(paste0(anc.bychr.dir,'/',AFN,'.A11.chr',chr,'.rds'))[filter.ids[fq[[2]]==chr],]
		anc.ea<-readRDS(paste0(anc.bychr.dir,'/',AFN,'.A12.chr',chr,'.rds'))[filter.ids[fq[[2]]==chr],]
		anc.sum<-anc.ee+anc.ea
		anc.sum[anc.sum<1]=1
		anc.ee<-anc.ee/anc.sum
		anc.ea<-anc.ea/anc.sum
		anc.aa<-1-anc.ee-anc.ea
		
		anc.anno<-imp.chrs[[chr]]
		anc.dis<-imp.anc.intv[[chr]]
		colnames(anc.dis)<-c('il','ir')
		anc.anno<-cbind(anc.anno,anc.dis)
		each.ind<-function(ind.id){
			anc.mat<-cbind(anc.ee[,ind.id],anc.ea[,ind.id],anc.aa[,ind.id])
			anc.stats<-apply(anc.mat,1,which.max)
			anc.bestc<-apply(anc.mat,1,max)
			anc.diff<-diff(anc.stats)
			anc.brkpts<-which(anc.diff!=0)
			np<-length(anc.brkpts)
			if(np==0){#no anc change
				anc.sp=1
				anc.ep=nrow(anc.mat)
				anc.se<-matrix(c(anc.anno[anc.sp,4],anc.anno[anc.ep,4]),nrow=1)
				anc.intvlen<-anc.se[,2]-anc.se[,1]
			}
			else{
				anc.sp<-c(1,anc.brkpts[1:np]+1)
				anc.ep<-c(anc.brkpts[1:np],nrow(anc.mat))
				anc.se<-cbind(anc.anno[anc.sp,4],anc.anno[anc.ep,4])
				anc.intvlen<-anc.se[,2]-anc.se[,1]
			}
			itv.num<-nrow(anc.se)
			mask.int<-c()
			for(intv.id in 1:itv.num){
				intv.st=anc.sp[intv.id]#save index start
				intv.ed=anc.ep[intv.id]#save index end
				wT<-anc.bestc[intv.st:intv.ed]<0.9
				if(any(wT)){
					ww<-which(wT)
					if(anc.intvlen[intv.id]<=dis.thres*2){
						anc.stats[intv.st:intv.ed]=NA
						mask.int<-rbind(mask.int,c(intv.st+ww[1]-1,intv.st,intv.ed))
					}
					else{
						for(j in ww){
							j.pos=intv.st+j-1
							anc.pos.l<-anc.anno[j.pos,'il']
							anc.pos.r<-anc.anno[j.pos,'ir']
							if(anc.pos.l>=intv.st){
								l.pos=anc.pos.l
							}
							else{
								l.pos=intv.st
							}
							if(anc.pos.r<=intv.ed){
								r.pos=anc.pos.r
							}
							else{
								r.pos=intv.ed
							}
							mask.int<-rbind(mask.int,c(j.pos,l.pos,r.pos))
						}
					}
				}
			}
			zz<-anc.stats
			maskind=rep(1,nrow(anc.mat))
			if(is.null(mask.int)){
				return(maskind)
			}
			if(nrow(mask.int)!=0){
				for(i in 1:nrow(mask.int)){
					maskind[mask.int[i,2]:mask.int[i,3]]=0
				}
			}
			return(maskind)	
		}
		ec<-do.call(cbind,mclapply(1:ncol(anc.ee),each.ind,mc.cores=mcc))
		colnames(ec)<-colnames(anc.ee)
		saveRDS(ec,paste0(outpdir,'/chr',chr,'_mask.rds'))
	}
	mask.chrs<-lapply(1:22,each_chr)
}

#functions below are used for ploting (inspecting the ancestry segments before/after running the masking)
#gen.bar.cp.mat.X<-function(chr){
#	imp.chr.anno<-imp.chrs[[chr]]
#	n.s<-nrow(imp.chr.anno)
#	imp.width<-diff(imp.chr.anno[,4])
#	imp.xb<-c(imp.chr.anno[1,4],imp.chr.anno[2:n.s,4]-imp.width/2)
#	imp.xt<-c(imp.chr.anno[1:(n.s-1),4]+imp.width/2,imp.chr.anno[n.s,4])
#	imp.x<-data.frame(bot=imp.xb,top=imp.xt)
#}
#
#if(!exists('imp.anno.chrs')){
#	imp.anno.chrs<-lapply(1:22,gen.bar.cp.mat.X)
#}
#gen.bar.cpss.check<-function(mask.info,shift=0){
#   imp.eu.yt<-mask.info$EE+shift
#   imp.uc.yt<-imp.eu.yt+mask.info$EA
#   heu<-cbind(imp.eu.yt,imp.uc.yt)
#	heu
#}
#
#get.chr.ind.check<-function(mask.list,chr,id,outdir='check_anc_filter'){
#   if(!dir.exists(outdir)){
#       dir.create(outdir)
#   }
#  	bar.dt<-gen.bar.cpss.check(mask.list$dt)
#	bar.dt1<-gen.bar.cpss.check(mask.list$dt,shift=1)
#	Xt<-imp.anno.chrs[[chr]]
#   png(paste0(outdir,'/chr',chr,'.',id,'.check.png'),width=4.875,height=2.2,units='in',res=1200,pointsize=4)
#	plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0('chr',chr),xlim=range(Xt[,1]),ylim=c(0,3.5),bty='n')
#
#	delt<-0.2
#   rect(Xt[,1],0+delt,Xt[,2],bar.dt[,1]+delt,col='red',border=NA)
#   rect(Xt[,1],bar.dt[,1]+delt,Xt[,2],bar.dt[,2]+delt,col='orange',border=NA)
#   rect(Xt[,1],bar.dt[,2]+delt,Xt[,2],1+delt,col='yellow',border=NA)
#
#	red.col<-rep('red',nrow(Xt))
#	org.col<-rep('orange',nrow(Xt))
#	yel.col<-rep('yellow',nrow(Xt))
#
#	for(i in 1:nrow(mask.list[[2]])){
#		red.col[mask.list[[2]][i,2]:mask.list[[2]][i,3]]='grey'
#		org.col[mask.list[[2]][i,2]:mask.list[[2]][i,3]]='grey'
#		yel.col[mask.list[[2]][i,2]:mask.list[[2]][i,3]]='grey'
#	}
#   rect(Xt[,1],1+2*delt,Xt[,2],bar.dt1[,1]+2*delt,col=red.col,border=NA)
#   rect(Xt[,1],bar.dt1[,1]+2*delt,Xt[,2],bar.dt1[,2]+2*delt,col=org.col,border=NA)
#   rect(Xt[,1],bar.dt1[,2]+2*delt,Xt[,2],2+2*delt,col=yel.col,border=NA)
#	axis(side=1)
#	legend('top',bty='n',ncol=4,fill=c('red','orange','yellow','grey'),legend=c('European','mixed','African','filtered'))
#   dev.off()
#}
#######################
out.mcgeno.dir<-paste0(out.geno.dir,'/mean_center_geno/')
if(!dir.exists(out.mcgeno.dir)){
	dir.create(out.mcgeno.dir)
}
mean.center.geno<-function(chr){
	#AFN<-basename(ADMIXINDFILE)
	geno.info<-readRDS(paste0(geno.bychr.dir,'/',AFN,'.rawgeno.chr',chr,'.rds'))
	anc.A11<-readRDS(paste0(anc.bychr.dir,'/',AFN,'.A11.chr',chr,'.rds'))
	anc.A12<-readRDS(paste0(anc.bychr.dir,'/',AFN,'.A12.chr',chr,'.rds'))
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
	#ANC2.global<-colMeans(ANC2)[filter.ids[fq[[2]]==chr]]#estimate global ancestry for POP2
	saveRDS(ANC2[filter.ids[fq[[2]]==chr],],paste0(out.mcgeno.dir,'/global_ANC_pop2.chr',chr,'.rds'))

	geno.P11<-as.matrix(readRDS(paste0(geno.bychr.dir,'/',AFN,'.pop11.chr',chr,'.rds')))
	geno.P12<-as.matrix(readRDS(paste0(geno.bychr.dir,'/',AFN,'.pop12.chr',chr,'.rds')))
	geno.P21<-as.matrix(readRDS(paste0(geno.bychr.dir,'/',AFN,'.pop21.rds',chr,'.rds')))
	geno.P22<-as.matrix(readRDS(paste0(geno.bychr.dir,'/',AFN,'.pop22.rds',chr,'.rds')))
		
	geno.P1<-geno.P11+geno.P12-(2*fq$pop1.freq*anc.A11+fq$pop1.freq*anc.A12)
	geno.P2<-geno.P22+geno.P21-(2*fq$pop2.freq*anc.A22+fq$pop2.freq*anc.A12)
	
	saveRDS(geno.P1[filter.ids[fq[[2]]==chr],],paste0(out.mcgeno.dir,'/',AFN,'.pop1.chr',chr,'.rds'))
	saveRDS(geno.P2[filter.ids[fq[[2]]==chr],],paste0(out.mcgeno.dir,'/',AFN,'.pop2.chr',chr,'.rds'))
}

run.meangeno<-function(){
	lapply(1:22,mean.center.geno)
	gl.anc<-do.call(rbind,lapply(1:22,function(chr)readRDS(paste0(out.mcgeno.dir,'/global_ANC_pop2.chr',chr,'.rds'))))
	glol.anc<-colMeans(gl.anc)
	saveRDS(glol.anc,paste0(out.geno.dir,'/all_inds_global_pop2_ancestry.rds'))
}

align.beta<-function(beta.df.all,chr){#beta.df.all should be a data frame
	beta.df<-beta.df.all[beta.df.all$chr==chr,]
	beta.ids.1<-paste0(beta.df$chr,':',beta.df$pos_grch37,':',beta.df$other_allele,':',beta.df$effect_allele)
	beta.ids.2<-paste0(beta.df$chr,':',beta.df$pos_grch37,':',beta.df$effect_allele,':',beta.df$other_allele)
	#fq.maf.snp<-paste0(fq.maf[[2]],':',fq.maf[[4]],':',fq.maf[[5]],':',fq.maf[[6]])
	fq.maf.snp<-paste0(imp.chrs[[chr]][[2]],':',imp.chrs[[chr]][[4]],':',imp.chrs[[chr]][[5]],':',imp.chrs[[chr]][[6]])
	
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
	each.chr<-function(chr){
		fbb<-align.beta(betas,chr)
		beta<-fbb$beta.ndf$BETA
		fq.ids<-fbb$fq.beta.ids
		g.p1<-readRDS(paste0(out.mcgeno.dir,'/',AFN,'.pop1.chr',chr,'.rds'))[fq.ids,]
		g.p2<-readRDS(paste0(out.mcgeno.dir,'/',AFN,'.pop2.chr',chr,'.rds'))[fq.ids,]

		#####masking the genomic regions to prevent those regions involved in PGS calculation
		f.chr<-readRDS(paste0(outpdir,'/chr',chr,'_mask.rds'))[fq.ids,]
		g.p1[f.chr==0]=0
		g.p2[f.chr==0]=0
		##############################
		pgs.p1<-t(g.p1)%*%matrix(beta,ncol=1)
		pgs.p2<-t(g.p2)%*%matrix(beta,ncol=1)
		pgs.chr<-data.frame(pgs.pop1=pgs.p1,pgs.pop2=pgs.p2)
		if(is.na(outpgs.prefix)){
			saveRDS(pgs,paste0(out.mcgeno.dir,'/',basename(betafile),'.admixpgs.chr',chr,'.rds'))
		}
		else{
			saveRDS(pgs,paste0(out.mcgeno.dir,'/',outpgs.prefix,'.chr',chr,'.rds'))
		}
		return(pgs.chr)
	}
	pgs.chrs<-lapply(1:22,each.chr)
	pgs.all<-Reduce('+',pgs.chrs)
	if(is.na(outpgs.prefix)){
		saveRDS(pgs.all,paste0(output.dir,'/',basename(betafile),'.admixpgs.rds'))
	}
	else{
		saveRDS(pgs.all,paste0(out.mcgeno.dir,'/',outpgs.prefix,'.rds'))
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
