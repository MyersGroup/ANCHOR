library(parallel)

if(!exists('imp.anno')){
	imp.anno<-readRDS('imp_anno_withmaf.rds')
	imp.order<-1:nrow(imp.anno)
	imp.chr<-mclapply(1:22,function(i)imp.order[!is.na(imp.anno$maf) &imp.anno[,2]==i & imp.anno$maf>=0.01 ],mc.cores=6)
}

if(!exists('imp.chrs')){
	imp.chrs<-lapply(1:22,function(i)imp.anno[imp.anno[,2]==i,])
	imp.chrs<-lapply(imp.chrs,function(x)x[!is.na(x$maf) & x$maf>=0.01,])
}

#This data object was obtained by applying the main function get.segment.wd
#imp.anc.5mb<-lapply(1:22,get.segment.wd,imp.chrs=imp.chrs,dis=2.5e6)
if(!exists('imp.anc.5mb')){
	imp.anc.5mb<-readRDS('imp_anc_5mb_dis.rds')
}

#dir.anc is the path where the local ancestry matrix (split by chromosomes) is stored
dir.anc<-'imp_anc'

#Af_anc is the global african ancestry for each individual
if(!exists('af.anc')){
	af.anc<-readRDS('African_ancestry_8003.rds')
}

get.segment.wd<-function(imp.chrs,chr,dis=2.5e6){
#get the 5MB window for all variants)
#We recommend to use 5MB (2.5MB towards each side from the centre)
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
	lz<-mclapply(1:nrow(g.chr),cpc,mc.cores=12)
	lz.dt<-do.call(rbind,lz)	
}

gen.bar.cp.mat.X<-function(chr){
    imp.chr.anno<-imp.anno[imp.chr[[chr]],]
    n.s<-nrow(imp.chr.anno)
    imp.width<-diff(imp.chr.anno[,4])
    imp.xb<-c(imp.chr.anno[1,4],imp.chr.anno[2:n.s,4]-imp.width/2)
    imp.xt<-c(imp.chr.anno[1:(n.s-1),4]+imp.width/2,imp.chr.anno[n.s,4])
	imp.x<-data.frame(bot=imp.xb,top=imp.xt)
}

#The object imp.anno.chrs was obtained by running
#   imp.anno.chrs<-lapply(1:22,gen.bar.cp.mat.X)

if(!exists('imp.anno.chrs')){
	imp.anno.chrs<-readRDS('var_x_position.rds')
}

#This function is used for generating the object for ploting the masked region
gen.bar.cpss<-function(chr,af.range=c(0.45,0.52)){
	bin.inds<-which(af.anc>af.range[1] & af.anc<=af.range[2])
	sf<-length(bin.inds)
	sf.q<-0:(sf-1)
	unphased.anc.copy<-readRDS(paste0(dir.anc,'/sample_all_8003_chr',chr,'.rds'))
	#### unphased.anc.copy is N by (A1A1+A1A2) local ancestry matrix  where A1A1 is the homogenous ancestry from pop1 and A1A2 is heterozygous ancestry,and A1A1+A1A2+A2A2=1, the local 
	unphased.anc.copy[unphased.anc.copy<0]=0
   	iac<-imp.anno.chrs[[chr]]
	n.s<-nrow(iac) 
	gen.bar.cp.mat<-function(id,shift=0){
    	imp.eu.yt<-unphased.anc.copy[,(id-1)*2+1]+shift
    	imp.uc.yt<-imp.eu.yt+unphased.anc.copy[,(id-1)*2+2]
    	heu<-cbind(imp.eu.yt,imp.uc.yt)
		#attr(heu,'shift')<-shift
		heu
	}
	gbc<-do.call(cbind,mclapply(sf.q,function(shift){
		g.id<-bin.inds[shift+1]
		cm<-gen.bar.cp.mat(g.id,shift)
	},mc.cores=12))
	ge<-gbc	
}
#
#get.chr.ind<-function(chr,outdir='anc_bin_plot2'){
#    if(!dir.exists(outdir)){
#        dir.create(outdir)
#    }
#   	bar.dt<-gen.bar.cpss(chr)
#	Xt<-imp.anno.chrs[[chr]]
#    png(paste0(outdir,'/chr',chr,'.gap.png'),width=4.875,height=27,units='in',res=1200,pointsize=4)
#	#plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0('chr',chr),xlim=range(Xt[,1]),ylim=c(0,ncol(bar.dt)/2+5),bty='n')
#	plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0('chr',chr),xlim=range(Xt[,1]),ylim=c(0,ncol(bar.dt)/2+10),bty='n')
#    tn<-ncol(bar.dt)/2
#	for(i in 1:tn){
#		delt<-(i-1)*0.1
#    	rect(Xt[,1],i-1+delt,Xt[,2],bar.dt[,(i-1)*2+1]+delt,col='red',border=NA)
#    	rect(Xt[,1],bar.dt[,(i-1)*2+1]+delt,Xt[,2],bar.dt[,(i-1)*2+2]+delt,col='orange',border=NA)
#    	rect(Xt[,1],bar.dt[,(i-1)*2+2]+delt,Xt[,2],i+delt,col='yellow',border=NA)
#	}
#    legend('top',bty='n',ncol=3,fill=c('red','orange','yellow'),legend=c('European','mixed','African'))
#    dev.off()
#}
##
#get.chr.ind.filter<-function(chr,af.range=c(0.45,0.52),outdir='anc_bin_plot_filter'){
#    if(!dir.exists(outdir)){
#        dir.create(outdir)
#    }
#	filter.dir<-'imp_mask_region_5mb'
#	filter.file<-paste0(filter.dir,'/chr',chr,'_5mb.rds')
#	filter.mat<-readRDS(filter.file)
#	bin.inds<-which(af.anc>af.range[1] & af.anc<=af.range[2])
#	filter.mat<-filter.mat[,bin.inds]
#   	bar.dt<-gen.bar.cpss(chr,af.range)
#	Xt<-imp.anno.chrs[[chr]]
#	xt.n<-nrow(Xt)
#    png(paste0(outdir,'/chr',chr,'.gap.png'),width=4.875,height=27,units='in',res=1200,pointsize=4)
#	#plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0('chr',chr),xlim=range(Xt[,1]),ylim=c(0,ncol(bar.dt)/2+5),bty='n')
#	plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0('chr',chr),xlim=range(Xt[,1]),ylim=c(0,ncol(bar.dt)/2+10),bty='n')
#    tn<-ncol(bar.dt)/2
#	for(i in 1:tn){
#		delt<-(i-1)*0.1
#
#		col1<-rep('red',xt.n)
#		col2<-rep('orange',xt.n)
#		col3<-rep('yellow',xt.n)
#		
#		col1[filter.mat[,i]==0]<-'grey'		
#		col2[filter.mat[,i]==0]<-'grey'		
#		col3[filter.mat[,i]==0]<-'grey'		
#
#    	rect(Xt[,1],i-1+delt,Xt[,2],bar.dt[,(i-1)*2+1]+delt,col=col1,border=NA)
#    	rect(Xt[,1],bar.dt[,(i-1)*2+1]+delt,Xt[,2],bar.dt[,(i-1)*2+2]+delt,col=col2,border=NA)
#    	rect(Xt[,1],bar.dt[,(i-1)*2+2]+delt,Xt[,2],i+delt,col=col3,border=NA)
#	}
#    legend('top',bty='n',ncol=4,fill=c('red','orange','yellow','grey'),legend=c('European','mixed','African','filter'))
#    dev.off()
#}

#run.all.ids<-function(){
#    mclapply(1:22,get.chr.ind.filter,mc.cores=12)
#}


anc.chr.filter.all<-function(chr,ind.id=NULL,dis.thres=2.5e6,mcc=12,outpdir='imp_mask_region_5mb/'){
	if(!dir.exists(outpdir)){
		dir.create(outpdir)
	}
	uac<-readRDS(paste0(dir.anc,'/sample_all_8003_chr',chr,'.rds'))
	anc.ee<-uac[,seq(1,ncol(uac),by=2)]
	anc.ea<-uac[,seq(2,ncol(uac),by=2)]
	anc.sum<-anc.ee+anc.ea
	anc.sum[anc.sum<1]=1
	anc.ee<-anc.ee/anc.sum
	anc.ea<-anc.ea/anc.sum
	anc.aa<-1-anc.ee-anc.ea
	anc.anno<-imp.chrs[[chr]]
	anc.dis5<-imp.anc.5mb[[chr]]
	anc.anno<-cbind(anc.anno,anc.dis5)
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
						anc.pos.l<-anc.anno[j.pos,9]
						anc.pos.r<-anc.anno[j.pos,10]
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
		#code below is used as input for function get.chr.ind.check 
		#rd<-data.frame(mask=anc.stats,bestcall=anc.bestc)
		#rd<-cbind(rd,anc.mat)
		#colnames(mask.int)<-c('cur.idx','left.idx','right.idx')
		#colnames(rd)<-c('Ancestry','best.call','EE','EA','AA')
		#rd.list<-list(dt=rd,mask_id=mask.int)
	}
#	#test chr22: 43:afanc:0.45, 37,afanc:0.34, 384, afanc:0.0.75
	if(!is.null(ind.id)){
		ec<-each.ind(ind.id)
	}
	else{#all ind
		ec<-do.call(cbind,mclapply(1:(ncol(uac)/2),each.ind,mc.cores=mcc))
		colnames(ec)<-names(af.anc)
		ec
		saveRDS(ec,paste0(outpdir,'/chr',chr,'_5mb.rds'))
	}
	return(ec)
}
#
#gen.bar.cpss.check<-function(mask.info,shift=0){
#    imp.eu.yt<-mask.info$EE+shift
#    imp.uc.yt<-imp.eu.yt+mask.info$EA
#    heu<-cbind(imp.eu.yt,imp.uc.yt)
#	heu
#}
#
#get.chr.ind.check<-function(mask.list,chr,id,outdir='check_anc_filter'){
#    if(!dir.exists(outdir)){
#        dir.create(outdir)
#    }
#   	bar.dt<-gen.bar.cpss.check(mask.list$dt)
#	bar.dt1<-gen.bar.cpss.check(mask.list$dt,shift=1)
#	Xt<-imp.anno.chrs[[chr]]
#    png(paste0(outdir,'/chr',chr,'.',id,'.check.png'),width=4.875,height=2.2,units='in',res=1200,pointsize=4)
#	plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0('chr',chr),xlim=range(Xt[,1]),ylim=c(0,3.5),bty='n')
#		
#	delt<-0.2
#    rect(Xt[,1],0+delt,Xt[,2],bar.dt[,1]+delt,col='red',border=NA)
#    rect(Xt[,1],bar.dt[,1]+delt,Xt[,2],bar.dt[,2]+delt,col='orange',border=NA)
#    rect(Xt[,1],bar.dt[,2]+delt,Xt[,2],1+delt,col='yellow',border=NA)
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
#    rect(Xt[,1],1+2*delt,Xt[,2],bar.dt1[,1]+2*delt,col=red.col,border=NA)
#    rect(Xt[,1],bar.dt1[,1]+2*delt,Xt[,2],bar.dt1[,2]+2*delt,col=org.col,border=NA)
#    rect(Xt[,1],bar.dt1[,2]+2*delt,Xt[,2],2+2*delt,col=yel.col,border=NA)
#	axis(side=1)
#	legend('top',bty='n',ncol=4,fill=c('red','orange','yellow','grey'),legend=c('European','mixed','African','filtered'))
#    dev.off()
#}
