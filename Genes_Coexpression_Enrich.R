########################################################  
				 # P2P correlation #
########################################################
colCors = function(x, y) {
	sqr = function(x) x*x
	if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
	stop("Please supply two matrices of equal size.")
	x   = sweep(x, 2, colMeans(x))
	y   = sweep(y, 2, colMeans(y))
	cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
	return(cor)
}
#######################################################  
				# Match Gene index # 
#######################################################
library(data.table)
sigfile = "/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt"
allfile = "/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt"
sig = read.table(sigfile,head=T,sep="\t",stringsAsFactors=F)
all = fread(allfile,data.table=F,head=T,sep="\t",stringsAsFactors=F)
promprb = read.table("/shares/compbio/Group-Yang/y.wu/interaction/data/promoter_mprobes.list",head=F)
ensem = read.delim("/shares/compbio/Group-Yang/y.wu/interaction/data/gencode_hg19.txt",head=F)
express = fread("/shares/compbio/Group-Yang/GTExV7/RNAseq/Tissue_Specific_Phen/GTExV7_gene_tpm_Whole_Blood.phen",head=T,data.table=F)[,-c(1:2)]
ensidx = match(colnames(express),ensem$V1,nomatch=0)
ensem = ensem[ensidx,]
##############################
		 # QC sig p2p 
##############################
pindx = match(sig$Outco_ID,promprb$V1,nomatch=0)
p2psig = sig[which(pindx!=0),]
indx = which(p2psig$Expo_Gene!=p2psig$Outco_Gene)
p2psig = p2psig[indx,]
sigdist = abs(p2psig$Expo_bp-p2psig$Outco_bp)
sigg1 = p2psig$Expo_Gene
sigg2 = p2psig$Outco_Gene
pairsig = paste(sigg1,sigg2,sep="-")
dupidx = which(duplicated(pairsig)==1)
sigg1 = sigg1[-dupidx]; sigg2 = sigg2[-dupidx]
sigdist = sigdist[-dupidx]
sigg1indx = match(sigg1,ensem$V9,nomatch=0)
sigg2indx = match(sigg2,ensem$V9,nomatch=0)
zindx = which(sigg1indx==0 | sigg2indx==0)
sigg1 = as.character(ensem$V1[sigg1indx[-zindx]])
sigg2 = as.character(ensem$V1[sigg2indx[-zindx]])
sigdist = sigdist[-zindx]
##############################
		 # QC all p2p
##############################
pindx = match(all$Outco_ID,promprb$V1,nomatch=0)
p2pall = all[which(pindx!=0),]
indx = which(p2pall$Expo_Gene!=p2pall$Outco_Gene)
p2pall = p2pall[indx,]
alldist = abs(p2pall$Expo_bp-p2pall$Outco_bp)
allg1 = p2pall$Expo_Gene
allg2 = p2pall$Outco_Gene
pairall = paste(allg1,allg2,sep="-")
dupidx = which(duplicated(pairall)==1)
allg1 = allg1[-dupidx]; allg2 = allg2[-dupidx]
alldist = alldist[-dupidx]
g1indx = match(allg1,ensem$V9,nomatch=0)
g2indx = match(allg2,ensem$V9,nomatch=0)
zindx = which(g1indx==0 | g2indx==0)
allg1 = as.character(ensem$V1[g1indx[-zindx]])
allg2 = as.character(ensem$V1[g2indx[-zindx]])
alldist = alldist[-zindx]
#######################################################  
				 # Matching distance
#######################################################
bin = 50
itvl=(max(alldist)-min(alldist))/bin
varbin = seq(itvl, max(alldist), by=itvl)
nvarbin = length(varbin)
varindx = ceiling(sigdist/itvl)
nperitvl=c()
for(i in 1:nvarbin)
{
    nperitvl=c(nperitvl,length(which(varindx==i)))
}

splsets<-matrix(list(), nrow=nvarbin, ncol=1)
for(i in 1:nvarbin)
{
    splsets[[i]]=which((alldist>(varbin[i]-itvl)) & (alldist<=varbin[i]))
}
#################### sampling #########################
cornm=c();
for(i in 1:1000){
	print(i); smplmprobes=c()
    for(j in 1:nvarbin)
    {
        smplmprobes=c(smplmprobes,sample(splsets[[j]],nperitvl[j],replace = FALSE))
    }
	cornull=colCors(as.matrix(express[,allg1[smplmprobes]]),as.matrix(express[,allg2[smplmprobes]]))
	cornm=c(cornm,mean(cornull,na.rm=T))
}
#######################################################  
					# Observed cor  
#######################################################
corobsv=colCors(as.matrix(express[,sigg1]),as.matrix(express[,sigg2]))

cornm=c();
for(i in 1:1000) {
	print(i)
	ranidx=sample(1:length(allg1),length(sigg1))
	cornull=colCors(as.matrix(express[,allg1[ranidx]]),as.matrix(express[,allg2[ranidx]]))
	cornm=c(cornm,mean(cornull,na.rm=T))
}
write.table(cornm,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/coexp/cor_null_1000.txt",row=F,col=F,sep="\t",qu=F)
####################################################    
				  # gene overlap     
####################################################   
PIR=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/eqtl/res/PIR_sig_heidi.txt",head=T,sep="\t")
m2m=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt",head=T,sep="\t")
m2g=read.table("/shares/compbio/Group-Yang/y.wu/omic2/result/m2gsmr_pass_heidi.smr",head=T)
pairPIR=paste(PIR$mprobe,PIR$gene,sep="-")
pairm2m=paste(m2m$Outco_ID,m2m$Expo_Gene,sep="-")
pairm2g=paste(m2g$Expo_ID,m2g$Outco_Gene,sep="-")
pairm2m=pairm2m[!duplicated(pairm2m)]
idx=match(pairPIR,pairm2g,nomatch=0)
sig=length(which(idx!=0))

nullall=c()
for(r in 1:1000){
	PIR=read.table(paste0("/shares/compbio/Group-Yang/y.wu/interaction/P2E/eqtl/jobs/samplePIR/PIR_null_",r,".txt"),head=T,sep="\t")
	pairPIR=paste(PIR$mprobe,PIR$gene,sep="-")
	idx=match(pairPIR,pairm2g,nomatch=0)
	null=length(which(idx!=0))
	nullall=c(nullall,null)
}
nullmean=mean(nullall)
nullsd=sd(nullall)

write.table(nullall,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/eqtl/res/null_smr_overlapgenes_1000.txt",row=F,col=F,sep="\t",qu=F)
t=qt(c(.025, .975), df=length(nullall)-1)
sd(nullall)*t/sqrt(500)
#######################################################  
					  # end #   
#######################################################