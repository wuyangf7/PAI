###		1. gene activity enrichment with all genes as controls
###		2. gene activity enrichment for Hi-C prediction from Jung et al.
####################################################################
	   # Expression Enrichment with all genes as controls #
####################################################################
library(data.table)
smr=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t",stringsAsFactors=F)
rpkm=fread("gunzip -c /shares/compbio/Group-Yang/GTExV7/RNAseq/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz",skip=2,head=T,data.table=F)
rpkm=rpkm[,c(1:2,55)]
nonindx=which(rpkm[,3]<=0.1)
##################  Matching  #####################
genesig=as.character(unique(smr$Expo_Gene))
gsigindx=match(genesig,rpkm[,2],nomatch=0)
genesig=genesig[which(gsigindx!=0)]
##################  Partition  #####################
inactive=rpkm[nonindx,2]; active=rpkm[-nonindx,]
active[,3]=log2(active[,3]); 
quantiles=summary(active[,3])
fst=active[which(active[,3]>quantiles[1] & active[,3]<quantiles[2]),2]
snd=active[which(active[,3]>quantiles[2] & active[,3]<quantiles[3]),2]
trd=active[which(active[,3]>quantiles[3] & active[,3]<quantiles[5]),2]
fth=active[which(active[,3]>quantiles[5] & active[,3]<quantiles[6]),2]
##################     Sig.    ######################
sigidxn=match(genesig,inactive,nomatch=0)
sigidx1=match(genesig,fst,nomatch=0)
sigidx2=match(genesig,snd,nomatch=0)
sigidx3=match(genesig,trd,nomatch=0)
sigidx4=match(genesig,fth,nomatch=0)
psign=length(which(sigidxn!=0))/length(genesig)
psig1=length(which(sigidx1!=0))/length(genesig)
psig2=length(which(sigidx2!=0))/length(genesig)
psig3=length(which(sigidx3!=0))/length(genesig)
psig4=length(which(sigidx4!=0))/length(genesig)
psigSMR = c(psign,psig1,psig2,psig3,psig4)
##################    left    #######################
null=rpkm[-gsigindx,]
pnull = c()
for(i in 1:1000){
	set.seed(i)
	sidx = sample(1:nrow(null),length(genesig))
	generdm = null[sidx,2]
	nonidxn=match(generdm,inactive,nomatch=0)
	nonidx1=match(generdm,fst,nomatch=0)
	nonidx2=match(generdm,snd,nomatch=0)
	nonidx3=match(generdm,trd,nomatch=0)
	nonidx4=match(generdm,fth,nomatch=0)
	pnonn=length(which(nonidxn!=0))/length(generdm)
	pnon1=length(which(nonidx1!=0))/length(generdm)
	pnon2=length(which(nonidx2!=0))/length(generdm)
	pnon3=length(which(nonidx3!=0))/length(generdm)
	pnon4=length(which(nonidx4!=0))/length(generdm)
	pnon = c(pnonn,pnon1,pnon2,pnon3,pnon4)
	pnull = rbind(pnull,pnon)
}
pnonSMR = apply(pnull,2,mean); pnonsdSMR = apply(pnull,2,sd)
####################################################################
	   # Expression Enrichment for Hi-C from Jung et al. #
####################################################################
###		 download the supplementary to get the gene list
library(data.table); library(stringr)
loopGen = read.table("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/sigGM_genelist.txt")
gsigindx=match(loopGen$V1,rpkm[,2],nomatch=0)
genesig = as.character(loopGen$V1[which(gsigindx!=0)])
##################     Sig.    ######################
sigidxn=match(genesig,inactive,nomatch=0)
sigidx1=match(genesig,fst,nomatch=0)
sigidx2=match(genesig,snd,nomatch=0)
sigidx3=match(genesig,trd,nomatch=0)
sigidx4=match(genesig,fth,nomatch=0)
psign=length(which(sigidxn!=0))/length(genesig)
psig1=length(which(sigidx1!=0))/length(genesig)
psig2=length(which(sigidx2!=0))/length(genesig)
psig3=length(which(sigidx3!=0))/length(genesig)
psig4=length(which(sigidx4!=0))/length(genesig)
psigloop = c(psign,psig1,psig2,psig3,psig4)
##################    left    #######################
null=rpkm[-gsigindx,]
pnull = c()
for(i in 1:1000){
	set.seed(i)
	sidx = sample(1:nrow(null),length(genesig))
	generdm = null[sidx,2]
	nonidxn=match(generdm,inactive,nomatch=0)
	nonidx1=match(generdm,fst,nomatch=0)
	nonidx2=match(generdm,snd,nomatch=0)
	nonidx3=match(generdm,trd,nomatch=0)
	nonidx4=match(generdm,fth,nomatch=0)
	pnonn=length(which(nonidxn!=0))/length(generdm)
	pnon1=length(which(nonidx1!=0))/length(generdm)
	pnon2=length(which(nonidx2!=0))/length(generdm)
	pnon3=length(which(nonidx3!=0))/length(generdm)
	pnon4=length(which(nonidx4!=0))/length(generdm)
	pnon = c(pnonn,pnon1,pnon2,pnon3,pnon4)
	pnull = rbind(pnull,pnon)
}
pnonloop = apply(pnull,2,mean); pnonsdloop = apply(pnull,2,sd)
###################################################
#					# plot (18.10.19)		      #	
###################################################
dataSMR=cbind(psigSMR,pnonSMR);dataSMRse=cbind(c(0,0,0,0,0),pnonsdSMR)
dataloop=cbind(psigloop,pnonloop);dataloopse=cbind(c(0,0,0,0,0),pnonsdloop)
pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/new/blood_gene_activity_additional_submit.pdf",height=8,width=15)
par(mar=c(5,5,3,2),mfrow=c(1,2))
plot<-barplot(t(dataSMR),col=c("maroon","navy"),axes=TRUE, axisnames=F,ylim=c(0,0.75),cex.lab=1.7,beside=T,cex.names=1.7,ylab="Proportion of genes",xlab="Expression level groups",cex.axis=1.7,main="Pm-PAI genes",cex.main=1.7)
legend("topright",c("Sig. Genes","Control Genes"),cex=1.7,fill=c("maroon","navy"),bty="n")
names=c("Inactive","1st","2nd","3rd","4th")
axis(1, at=c(2, 5, 8, 11, 14), labels=names, srt=45,las=1, cex.axis=1.7)
segments(plot,t(dataSMR)+t(dataSMRse),plot,t(dataSMR)-t(dataSMRse), lwd=1.5)
arrows(plot,t(dataSMR)+t(dataSMRse),plot,t(dataSMR)-t(dataSMRse), angle=90, code=3, length=0.03)
mtext("a",side=3,line=0.5,adj=-0.15,cex=3.5,col="black",font=1)
###				plot loops enrichment
plot<-barplot(t(dataloop),col=c("maroon","navy"),axes=TRUE, axisnames=F,ylim=c(0,0.75),cex.lab=1.7,beside=T,cex.names=1.7,ylab="Proportion of genes",xlab="Expression level groups",cex.axis=1.7,main="Predicted target genes from PCHi-C",cex.main=1.7)
legend("topright",c("Sig. Genes","Control Genes"),cex=1.7,fill=c("maroon","navy"),bty="n")
names=c("Inactive","1st","2nd","3rd","4th")
axis(1, at=c(2, 5, 8, 11, 14), labels=names, srt=45,las=1, cex.axis=1.7)
segments(plot,t(dataloop)+t(dataloopse),plot,t(dataloop)-t(dataloopse), lwd=1.5)
arrows(plot,t(dataloop)+t(dataloopse),plot,t(dataloop)-t(dataloopse), angle=90, code=3, length=0.03)
mtext("b",side=3,line=0.5,adj=-0.15,cex=3.5,col="black",font=1)
dev.off()