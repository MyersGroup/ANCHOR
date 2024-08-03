library(parallel)

cc<-commandArgs(T)

if(length(cc)==0){
	print("We assume you would run this script within an R session")
	print("Please provide the following files:")
	print("ADMIXINDFIE=")
	print("OUTDIR=")
	print("HAPMIX_DATADIR=")
	print("ADMIXPOP")
	print("HAPMIX_MODE")
	print("output.dir")
	print("mcc=")
}
else{
	ADMIXINDFIE=cc[1]
	OUTDIR=cc[2]
	HAPMIX_DATADIR=cc[3]
	ADMIXPOP=cc[4]
	HAPMIX_MODE=cc[5]
	output.dir=cc[6]
	mcc=as.integer(cc[7])
}

if(!exists('mcc') || is.na(mcc)){
	mcc=12
}

#ADMIXINDFIE: individuals sample ID file
#OUTDIR: directory of hapmix output 
#ADMIXPOP:hapmix output file prefix
#HAPMIX_MODE: diploid as default 
#output.dir: store the output of this function 

column.map<-function(p1,h1,p2,h2){
    #p:{0,1} h:{0,1}
    col<-sum(c(p1,h1,p2,h2)*2^(3:0))+1
    return(col)
}

w<-function(p1,p2,h1,h2) column.map(p1,h1,p2,h2)

#pop1.geno<-function(prob16){#Don't use this function before mean-center genotype
#    pop1.1.col<-c(w(0,0,1,0),w(0,1,1,0),w(0,1,1,1),w(0,0,0,1),w(1,0,0,1),w(1,0,1,1))
#    pop1.2.col<-w(0,0,1,1)
#    pop1<-rowSums(prob16[,pop1.1.col])+2*prob16[,pop1.2.col]
#}
#
#pop2.geno<-function(prob16){#Don't use this function before mean-center genotype
#    pop2.1.col<-c(w(1,0,1,0),w(1,1,1,0),w(1,0,1,1),w(0,1,0,1),w(1,1,0,1),w(0,1,1,1))
#    pop2.2.col<-w(1,1,1,1)
#    pop2<-rowSums(prob16[,pop2.1.col])+2*prob16[,pop2.2.col]
#}

pop12.geno<-function(prob16){
	pop11.col<-c(w(0,0,1,0),w(0,0,0,1),w(0,0,1,1),w(0,0,1,1))
	pop12.col<-c(w(0,1,1,0),w(0,1,1,1),w(1,0,0,1),w(1,0,1,1))
	pop11.g<-rowSums(prob16[,pop11.col])
	pop12.g<-rowSums(prob16[,pop12.col])
	pop.g12<-data.frame(P11=pop11.g,P12=pop12.g)
}

pop21.geno<-function(prob16){
	pop22.col<-c(w(1,1,1,0),w(1,1,0,1),w(1,1,1,1),w(1,1,1,1))
	pop21.col<-c(w(1,0,1,0),w(1,0,1,1),w(0,1,0,1),w(0,1,1,1))
	pop22.g<-rowSums(prob16[,pop22.col])
	pop21.g<-rowSums(prob16[,pop21.col])
	pop.g21<-data.frame(P22=pop22.g,P21=pop21.g)
}

anc.prob<-function(prob16){
	prob.0.cols<-c(w(0,0,0,0),w(0,0,0,1),w(0,0,1,0),w(0,0,1,1))
	prob.2.cols<-c(w(1,1,0,0),w(1,1,0,1),w(1,1,1,0),w(1,1,1,1))
	prob.0<-rowSums(prob16[,prob.0.cols])
	prob.2<-rowSums(prob16[,prob.2.cols])
	prob.1<-1-prob.0-prob.2
	prob.1[prob.1<0]<-0
	pz<-data.frame(A11=prob.0,A12=prob.1,A22=prob.2)
	return(pz)
}

#read.imp<-function(ADMIXINDFILE,OUTDIR,ADMIXPOP,HAPMIX_MODE,output.dir,mcc=16){#don't use this function without mean-centering
#	#input arguments are exactly the same as input configuration file for Hapmix (.par file) except the output.dir, which specified by user 
#	if(!dir.exists(output.dir)){
#		dir.create(output.dir)
#	}
#	sample.ind<-read.table(ADMIXINDFILE,as.is=T)
#	n.sample<-nrow(sample.ind)-1
# 	if(!dir.exists(OUTDIR)){
#		stop(paste0('error: directory:',OUTDIR,' does not exist!'))
#	}
#    get.batch<-function(chr){
#        get.anc.geno<-do.call(cbind,lapply(0:n.sample,function(i){
#            fd<-read.table(paste0(OUTDIR,'/',ADMIXPOP,'.',HAPMIX_MODE,'.',i,'.',chr),as.is=T,header=F)
#            fd.p1<-pop1.geno(fd)
#            fd.p2<-pop2.geno(fd)
#            gn<-cbind(fd.p1,fd.p2)
#        }))
#    }
#    bgb<-do.call(rbind,mclapply(1:22,get.batch,mc.cores=mcc))
#    sp.ids<-sample.ind[,1]
#	colnames(bgb)<-paste0(c('pop1-','pop2-'),rep(sp.ids,each=2))
#    saveRDS(bgb,paste0(output.dir,'/',basename(ADMIXINDFILE),'.pop1_pop2.rds'))
#}

#This function use the input snp annotation file for hapmix to generate the variants annotation file

read.imp.anno<-function(HAPMIX_DATADIR,mcc=10){
    read.imp.chr.anno<-function(chr){
        chr.anno.org<-read.table(paste0(HAPMIX_DATADIR,'/chr',chr,'_snpfile.ORIG'),as.is=T,header=F)
        chr.anno<-read.table(paste0(HAPMIX_DATADIR,'/chr',chr,'_snpfile'),as.is=T,header=F)
		c1<-paste0(chr.anno.org[,2],':',chr.anno.org[,4],':',chr.anno.org[,5],':',chr.anno.org[,6])
		c2<-paste0(chr.anno[,2],':',chr.anno[,4],':',chr.anno[,5],':',chr.anno[,6])
		chr.anno.org$goodsnp<-c1%in%c2
		chr.anno.org
    }
    all.anno<-do.call(rbind,mclapply(1:22,read.imp.chr.anno,mc.cores=mcc))
}

#The function below use the same variable defined in the configuration file of hapmix:
#ADMIXINDFIE: individuals sample ID file
#OUTDIR: directory of hapmix output 
#ADMIXPOP:hapmix output file prefix
#HAPMIX_MODE: diploid as default 
#output.dir: store the output of this function 
 
#read.imp4<-function(ADMIXINDFILE,OUTDIR,ADMIXPOP,HAPMIX_DATADIR,HAPMIX_MODE,output.dir,mcc=16){
read.imp4<-function(){
	#input arguments are exactly the same as input configuration file for Hapmix (.par file) except the output.dir, which specified by user 
	if(!dir.exists(output.dir)){
		dir.create(output.dir)
	}
	sample.ind<-read.table(ADMIXINDFILE,as.is=T)
	n.sample<-nrow(sample.ind)-1
 	if(!dir.exists(OUTDIR)){
		stop(paste0('error: directory:',OUTDIR,' does not exist!'))
	}
   	geno.bychr.dir<-paste0(output.dir,'/imp_geno/')
	anc.bychr.dir<-paste0(output.dir,'/imp_anc/')
	if(!dir.exists(geno.bychr.dir)){
		dir.create(geno.bychr.dir)
	}
	if(!dir.exists(anc.bychr.dir)){
		dir.create(anc.bychr.dir)
	}
    sp.ids<-sample.ind[,1]
	ria<-read.imp.anno(HAPMIX_DATADIR,mcc)
	saveRDS(ria[ria$goodsnp,],paste0(output.dir,'/SNP_annot.rds'))
	
	get.batch<-function(chr){
        get.anc.geno<-do.call(cbind,lapply(0:n.sample,function(i){
            fd<-read.table(paste0(OUTDIR,'/',ADMIXPOP,'.',HAPMIX_MODE,'.',i,'.',chr),as.is=T,header=F)
            fd.pop12<-pop12.geno(fd)
            fd.pop21<-pop21.geno(fd)
			raw.geno<-rowSums(fd.pop12)+rowSums(fd.pop21)
			fd.anc<-anc.prob(fd)
            gn<-cbind(fd.pop12,fd.pop21,data.frame(geno=raw.geno),fd.anc)
        }))
		chr.anno<-ria[ria[[2]]==chr,]
		chr.ids<-chr.anno$goodsnp
		pp.names<-c('pop11','pop12','pop22','pop21')
		for(pp in 1:4){
			pp.sp<-get.anc.geno[chr.ids,8*(0:n.sample)+pp]
			colnames(pp.sp)<-sp.ids
    		saveRDS(pp.sp,paste0(geno.bychr.dir,'/',basename(ADMIXINDFILE),'.',pp.names[pp],'.chr',chr,'.rds'))
		}
		rg<-get.anc.geno[chr.ids,8*(0:n.sample)+5]
		colnames(rg)<-sp.ids
		saveRDS(rg,paste0(geno.bychr.dir,'/',basename(ADMIXINDFILE),'.rawgeno.chr',chr,'.rds'))
		aa.names<-c('A11','A12','A22')
		for(aa in 6:8){
			aa.sp<-get.anc.geno[chr.ids,8*(0:n.sample)+aa]
			colnames(aa.sp)<-sp.ids
			saveRDS(aa.sp,paste0(anc.bychr.dir,'/',basename(ADMIXINDFILE),'.',aa.names[aa-5],'.chr',chr,'.rds'))
		}
    }
    #bgb<-do.call(rbind,mclapply(1:22,get.batch,mc.cores=mcc))
	#bgb<-bgb[ria$goodsnp,]
	bgb<-mclapply(1:22,get.batch,mc.cores=mcc)
}

if(length(cc)>0){
	r4<-read.imp4()
	q(save='no')
}
