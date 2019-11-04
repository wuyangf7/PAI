###################################################################################
#####		show the Hi-C loops from Grubert et al. with PAIs             #####
#####		     show the TAD annotation from Rao et al.           	      #####
###################################################################################
#####################################################################
#############		Fig. 2d, MAD1L1 locus		#############
#####################################################################
#####	 Extract regional Hi-C data
file="/shares/compbio/Group-Yang/y.wu/interaction/data/HiC/HiC.GM12878.correlations.txt.gz"
zcat $file | awk '{if ($1==7 && $6>=1000000 && $6<=3000000 && $7>=1000000 && $7<=3000000) print $0}' > /shares/compbio/Group-Yang/y.wu/interaction/HiC/plotdata/MAD1L1.txt
#####	 Hi-C correlation plot
library(data.table)
tad = fread("gunzip -c /shares/compbio/Group-Yang/y.wu/interaction/data/HiC/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz",head=T)
map = fread("gunzip -c /shares/compbio/Group-Yang/y.wu/interaction/data/HiC/fragmentStartEnd.txt.gz",head=F)
value = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/plotdata/MAD1L1.txt",head=F)
chr = 7; start = 1000000; end = 3000000
indx=which(map$V1==chr & map$V4>=start & map$V4<=end)
dat=matrix(0,length(indx),length(indx))
colnames(dat)=rownames(dat)=map$V4[indx]
cindx=match(value$V6,colnames(dat),nomatch=0)
rindx=match(value$V7,rownames(dat),nomatch=0)
for(i in 1:dim(value)[1]){
	dat[cindx[i],rindx[i]]=value$V5[i]
	dat[rindx[i],cindx[i]]=value$V5[i]
}
dat[which(dat<=0.4)]=0
gap=(end-start)/4
idx1=which.min(abs(as.numeric(colnames(dat))-(start)))
idx2=which.min(abs(as.numeric(colnames(dat))-(start+gap)))
idx3=which.min(abs(as.numeric(colnames(dat))-(start+2*gap)))
idx4=which.min(abs(as.numeric(colnames(dat))-(start+3*gap)))
idx5=which.min(abs(as.numeric(colnames(dat))-end))
#####	  plot pdf file
pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/HiC_interaction_fig2d.pdf")
library(corrplot)
num = dim(dat)[1]
a<-corrplot(dat,is.corr=F,method="color",tl.pos="n",
cl.lim = NULL,cl.length = NULL,cl.cex = 0.8,cl.align.text="c",cl.offset=0.5,mar=c(0,4,4,2))
segments(0, 0, num, 0, col= 'black')
labelA=c(round(end/1e6,1),round((start+3*gap)/1e6,1),round((start+2*gap)/1e6,1),round((start+1*gap)/1e6,1),round(start/1e6,1))
labelB=c(round(start/1e6,1),round((start+gap)/1e6,1),round((start+2*gap)/1e6,1),round((start+3*gap)/1e6,1),round(end/1e6,1))
axis(side=2,at=c(num-idx1+1,num-idx2+1,num-idx3+1,num-idx4+1,num-idx5),labels=labelB,line=-0.2)
axis(side=3,at=c(idx5,idx4-1,idx3-1,idx2-1,idx1-1),labels=labelA,line=-3.15)
axis(side=4,at=c(idx5,idx4-1,idx3-1,idx2-1,idx1-1),labels=c("","","","",""),line=-3.65)
#####	 SMR predicted PAIs plot
smr=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t",stringsAsFactors=F)
idx=which(smr$Expo_Chr==chr & smr$Expo_bp>=start & smr$Expo_bp<=end & smr$Outco_Chr==chr & smr$Outco_bp>=start & smr$Outco_bp<=end)
pos1=as.numeric(smr$Expo_bp[idx])
pos2=as.numeric(smr$Outco_bp[idx])
for(i in 1:length(pos1)){
	output="*"
	idx1=which.min(abs(as.numeric(rownames(dat))-pos1[i]))
	idx2=which.min(abs(as.numeric(rownames(dat))-pos2[i]))
	text(idx2,num-idx1+1,output,cex=1.4,col="red")
}
#####	 TAD from Rao et al. 
tadidx = which(tad$chr1==chr & tad$x1>=start & tad$x2<=end)
tadpos1 = as.numeric(tad$x1[tadidx])
tadpos2 = as.numeric(tad$x2[tadidx])
for(i in 1:length(tadpos1)){
	idx1=which.min(abs(as.numeric(rownames(dat))-tadpos1[i]))
	idx2=which.min(abs(as.numeric(rownames(dat))-tadpos2[i]))
	rect(idx1,num-idx1+1,idx2,num-idx2+1,col=NA,border=T,lwd=2)	
}
#####################################################################
#############		Fig. S3a, RNASET2 locus		#############
#####################################################################
file="/shares/compbio/Group-Yang/y.wu/interaction/data/HiC/HiC.GM12878.correlations.txt.gz"
zcat $file | awk '{if ($1==6 && $6>=166740000 && $6<=168120000 && $7>=166740000 && $7<=168120000) print $0}' > /shares/compbio/Group-Yang/y.wu/interaction/HiC/plotdata/RNASET2.txt
#####	 Hi-C correlation plot
library(data.table)
tad = fread("gunzip -c /shares/compbio/Group-Yang/y.wu/interaction/data/HiC/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz",head=T)
map=fread("gunzip -c /shares/compbio/Group-Yang/y.wu/interaction/data/HiC/fragmentStartEnd.txt.gz",head=F)
value=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/HiC/RNASET2.txt",head=F)

chr=6;start=166740000;end=168120000
indx=which(map$V1==chr & map$V4>=start & map$V4<=end)
dat=matrix(0,length(indx),length(indx))
colnames(dat)=rownames(dat)=map$V4[indx]
cindx=match(value$V6,colnames(dat),nomatch=0)
rindx=match(value$V7,rownames(dat),nomatch=0)
for(i in 1:dim(value)[1]){
	dat[cindx[i],rindx[i]]=value$V5[i]
	dat[rindx[i],cindx[i]]=value$V5[i]
}
dat[which(dat<=0.4)]=0
gap=(end-start)/4
idx1=which.min(abs(as.numeric(colnames(dat))-(start)))
idx2=which.min(abs(as.numeric(colnames(dat))-(start+gap)))
idx3=which.min(abs(as.numeric(colnames(dat))-(start+2*gap)))
idx4=which.min(abs(as.numeric(colnames(dat))-(start+3*gap)))
idx5=which.min(abs(as.numeric(colnames(dat))-end))
#####	  plot pdf file
pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/HiC_interaction_figS3a.pdf")
library(corrplot)
num = dim(dat)[1]
a<-corrplot(dat,is.corr=F,method="color",tl.pos="n",
cl.lim = NULL,cl.length = NULL,cl.cex = 0.8,cl.align.text="c",cl.offset=0.5,mar=c(0,4,4,2))
segments(0, 0, num, 0, col= 'black')
labelA=c(round(end/1e6,1),round((start+3*gap)/1e6,1),round((start+2*gap)/1e6,1),round((start+1*gap)/1e6,1),round(start/1e6,1))
labelB=c(round(start/1e6,1),round((start+gap)/1e6,1),round((start+2*gap)/1e6,1),round((start+3*gap)/1e6,1),round(end/1e6,1))
axis(side=2,at=c(num-idx1+1,num-idx2+1,num-idx3+1,num-idx4+1,num-idx5),labels=labelB,line=-0.2)
axis(side=3,at=c(idx5,idx4-1,idx3-1,idx2-1,idx1-1),labels=labelA,line=-3.15)
axis(side=4,at=c(idx5,idx4-1,idx3-1,idx2-1,idx1-1),labels=c("","","","",""),line=-3.65)
#####	 SMR predicted PAIs plot
smr=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t",stringsAsFactors=F)
idx=which(smr$Expo_Chr==chr & smr$Expo_bp>=start & smr$Expo_bp<=end & smr$Outco_Chr==chr & smr$Outco_bp>=start & smr$Outco_bp<=end)
pos1=as.numeric(smr$Expo_bp[idx])
pos2=as.numeric(smr$Outco_bp[idx])
for(i in 1:length(pos1)){
	output="*"
	idx1=which.min(abs(as.numeric(rownames(dat))-pos1[i]))
	idx2=which.min(abs(as.numeric(rownames(dat))-pos2[i]))
	text(idx2,num-idx1+1,output,cex=1.4,col="red")
}
#####	 TAD from Rao et al. 
tadidx = which(tad$chr1==chr & tad$x1>=start & tad$x2<=end)
tadpos1 = as.numeric(tad$x1[tadidx])
tadpos2 = as.numeric(tad$x2[tadidx])
for(i in 1:length(tadpos1)){
	idx1=which.min(abs(as.numeric(rownames(dat))-tadpos1[i]))
	idx2=which.min(abs(as.numeric(rownames(dat))-tadpos2[i]))
	rect(idx1,num-idx1+1,idx2,num-idx2+1,col=NA,border=T,lwd=2)	
}
#####################################################################
#############		Fig. S3b, ABCB9 locus		#############
#####################################################################
file="/shares/compbio/Group-Yang/y.wu/interaction/data/HiC/HiC.GM12878.correlations.txt.gz"
zcat $file | awk '{if ($1==12 && $6>=123320000 && $6<=124130000 && $7>=123320000 && $7<=124130000) print $0}' > /shares/compbio/Group-Yang/y.wu/interaction/HiC/plotdata/ABCB9.txt
#####	 Hi-C correlation plot
library(data.table)
tad = fread("gunzip -c /shares/compbio/Group-Yang/y.wu/interaction/data/HiC/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz",head=T)
map=fread("gunzip -c /shares/compbio/Group-Yang/y.wu/interaction/data/HiC/fragmentStartEnd.txt.gz",head=F)
value=fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/plotdata/ABCB9.txt",head=F)
chr = 12;start = 123326611;end = 124128451
indx=which(map$V1==chr & map$V4>=start & map$V4<=end)
dat=matrix(0,length(indx),length(indx))
colnames(dat)=rownames(dat)=map$V4[indx]
cindx=match(value$V6,colnames(dat),nomatch=0)
rindx=match(value$V7,rownames(dat),nomatch=0)
for(i in 1:dim(value)[1]){
	dat[cindx[i],rindx[i]]=value$V5[i]
	dat[rindx[i],cindx[i]]=value$V5[i]
}
dat[which(dat<=0.4)]=0
gap=(end-start)/4
idx1=which.min(abs(as.numeric(colnames(dat))-(start)))
idx2=which.min(abs(as.numeric(colnames(dat))-(start+gap)))
idx3=which.min(abs(as.numeric(colnames(dat))-(start+2*gap)))
idx4=which.min(abs(as.numeric(colnames(dat))-(start+3*gap)))
idx5=which.min(abs(as.numeric(colnames(dat))-end))
#####	  plot pdf file
pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/HiC_interaction_figS3b.pdf")
library(corrplot)
num = dim(dat)[1]
a<-corrplot(dat,is.corr=F,method="color",tl.pos="n",
cl.lim = NULL,cl.length = NULL,cl.cex = 0.8,cl.align.text="c",cl.offset=0.5,mar=c(0,4,4,2))
segments(0, 0, num, 0, col= 'black')
labelA=c(round(end/1e6,1),round((start+3*gap)/1e6,1),round((start+2*gap)/1e6,1),round((start+1*gap)/1e6,1),round(start/1e6,1))
labelB=c(round(start/1e6,1),round((start+gap)/1e6,1),round((start+2*gap)/1e6,1),round((start+3*gap)/1e6,1),round(end/1e6,1))
axis(side=2,at=c(num-idx1+1,num-idx2+1,num-idx3+1,num-idx4+1,num-idx5),labels=labelB,line=-0.2)
axis(side=3,at=c(idx5,idx4-1,idx3-1,idx2-1,idx1-1),labels=labelA,line=-3.15)
axis(side=4,at=c(idx5,idx4-1,idx3-1,idx2-1,idx1-1),labels=c("","","","",""),line=-3.65)
#####	 SMR predicted PAIs plot
smr=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t",stringsAsFactors=F)
idx=which(smr$Expo_Chr==chr & smr$Expo_bp>=start & smr$Expo_bp<=end & smr$Outco_Chr==chr & smr$Outco_bp>=start & smr$Outco_bp<=end)
pos1=as.numeric(smr$Expo_bp[idx])
pos2=as.numeric(smr$Outco_bp[idx])
for(i in 1:length(pos1)){
	output="*"
	idx1=which.min(abs(as.numeric(rownames(dat))-pos1[i]))
	idx2=which.min(abs(as.numeric(rownames(dat))-pos2[i]))
	text(idx2,num-idx1+1,output,cex=1.4,col="red")
}
#####	 TAD from Rao et al. 
tadidx = which(tad$chr1==chr & tad$x1>=start & tad$x2<=end)
tadpos1 = as.numeric(tad$x1[tadidx])
tadpos2 = as.numeric(tad$x2[tadidx])
for(i in 1:length(tadpos1)){
	idx1=which.min(abs(as.numeric(rownames(dat))-tadpos1[i]))
	idx2=which.min(abs(as.numeric(rownames(dat))-tadpos2[i]))
	rect(idx1,num-idx1+1,idx2,num-idx2+1,col=NA,border=T,lwd=2)	
}
#####################################################################
#############			end			#############
#####################################################################
