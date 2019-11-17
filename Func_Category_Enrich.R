############################################################################
#       Compute the enrichment of PIDSs in 14 functional categories        # 
#                       from the Roadmap project                           #  
############################################################################
#########################################################################################
#                         1. Estimate the fold enrichment                               #  
#########################################################################################
sigfile = "/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt"
allfile = "/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt"
sigDNAm = unique(read.table(sigfile,header=T,sep="\t")[,c(6,5,8)])
allDNAm = unique(fread(allfile,head=T,data.table=F)[,c(6,5,8)])
mevar=read.table("/shares/compbio/Group-Yang/y.wu/interaction/data/lbc.vp.txt",head=F)
indx=match(allDNAm$Outco_ID,mevar$V1,nomatch=0)
allDNAm$var=mevar$V2[indx]
options(scipen=999)
sigindx=match(sigDNAm$Outco_ID,allDNAm$Outco_ID,nomatch=0)
sigvar=allDNAm$var[sigindx]
######          Partition distribution
bin=100
itvl=(max(allDNAm$var)-min(allDNAm$var))/bin
varbin = seq(itvl, max(allDNAm$var), by=itvl)
nvarbin = length(varbin)
varindx = ceiling(sigvar/itvl)
nperitvl=c()
for(i in 1:nvarbin)
{
    nperitvl=c(nperitvl,length(which(varindx==i)))
}
######          Distribution intervals
splsets<-matrix(list(), nrow=nvarbin, ncol=1)
for(i in 1:nvarbin)
{
    splsets[[i]]=which((allDNAm$var>(varbin[i]-itvl)) & (allDNAm$var<=varbin[i]))
}
######          Null PIDSs sampling
for(i in 1:1000){
    print(i)
    smplmprobes=c()
    for(j in 1:nvarbin)
    {
        set.seed(i)
        smplmprobes=c(smplmprobes,sample(splsets[[j]],nperitvl[j],replace = FALSE))
    }
    outfile=paste0("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/samples/null_mprobes_",i,".bim")
    write.table(allDNAm[smplmprobes,c(1:3)],outfile,row=F,col=F,sep="\t",quote=F)
} 
######          observed PIDSs formatting
write.table(sigDNAm,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.bim",row=F,col=F,sep="\t",quote=F)
################################################ 
################## annotation ##################
################################################
anno="/shares/compbio/Group-Yang/y.wu/omics/enrichment/script/anno14imb"
bim="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/samples/null_mprobes_"{TASK_ID}".bim"
out="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/jobs/null_mprobes_"{TASK_ID}
log="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/log/null_mprobes_"{TASK_ID}".log"
cmd="$anno --cat 1,2,3,4,5,6,7,8,9,10,11,12,13,14 --bim $bim --out $out >$log 2>&1"
qsubshcom "$cmd" 1 3G anno 2:00:00 "-J 1-1000"

######     sig. PIDSs 
anno="/shares/compbio/Group-Yang/y.wu/omics/enrichment/script/anno14imb"
$anno --cat 1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
--bim /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.bim \
--out /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr \
> /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.log 2>&1

######     All PIDSs
anno="/shares/compbio/Group-Yang/y.wu/omics/enrichment/script/anno14imb"
bim="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/allprobes.bim"
out="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/allprobes"
log="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/allprobes.log"
rsub "$anno --cat 1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
--bim $bim --out $out >$log 2>&1" @anno %20 T1 W1
########################################################
#            Compute the fold enrichment               #
######################################################## 
library(data.table)
simu=1000
expall=matrix(0,14,simu)
for(i in 1:simu){
    print(i)
    filename=paste0("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/jobs/null_mprobes_",i,".anno.txt")
    null=fread(filename,head=T)
    rowsum=rowSums(null[,2:dim(null)[2]])
    expt=rowsum/sum(rowsum)
    expall[,i]=expt
}
write.table(expall,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/null_mprobes.txt",row=F,col=F,sep="\t",qu=F)
#########################################################################################
#                                 2. Summarize and plot                                 #  
######################################################################################### 
alter=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.anno.txt",head=T)
siganno=rowSums(alter[,2:dim(alter)[2]])/sum(rowSums(alter[,2:dim(alter)[2]]))
expanno=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/null_mprobes.txt",head=F)
enrich=siganno/apply(expanno,1,mean)
p_value=c()
for(i in 1:14){
    p=length(which(expanno[i,]>siganno[i]))/dim(expanno)[2]
    p_value=c(p_value,p)
}
output = data.frame(alter[,1],enrich,p_value)
output[order(output$enrich),]
enrichse = apply(siganno/expanno,1,sd)
######      plot 
library(stringr)
colfile=unique(read.delim("/gpfs/gpfs01/polaris/Q0286/uqywu16/chromstate/annotation_25_imputed12marks.txt",head=T,sep="\t")$COLOR.CODE)
col=str_split_fixed(colfile,",",3)

pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/new/sigmprobes_enrich_smr_submit_revise.pdf",height=8,width=15)
par(mar=c(5,5,3,2),mfrow=c(1,2))

library(data.table)
obsvfile=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.anno.txt",head=T)
obsv=rowSums(obsvfile[,2:dim(obsvfile)[2]])/sum(rowSums(obsvfile[,2:dim(obsvfile)[2]]))
expfile=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/allprobes.anno.txt",head=T,data.table=F)
exp=rowSums(expfile[,2:dim(expfile)[2]])/sum(rowSums(expfile[,2:dim(expfile)[2]]))

data=cbind(obsv,exp)
names=c("TssA","Prom","Tx","TxWk","TxEn","EnhA","EnhW","DNase","ZNF/Rpts","Het","PromP","PromBiv","ReprPC","Quies")
plot<-barplot(t(data),col=c("cornflowerblue","forestgreen"),axes=TRUE, axisnames=F,ylim=c(0,0.35),cex.lab=1.7,cex.axis=1.7,beside=T,cex.names=1.7,ylab="Proportion of probes") 
legend("topright",c("PIDSs","All DNAm"),cex=1.7,fill=c("cornflowerblue","forestgreen"),bty="n")
text(apply(plot,2,mean),par("usr")[3]-0.005,labels=names,srt=45,adj=1,xpd=TRUE,cex=1.7)

names=c("TssA","Prom","Tx","TxWk","TxEn","EnhA","EnhW","DNase","ZNF/Rpts","Het","PromP","PromBiv","ReprPC","Quies")   
plot<-barplot(enrich,col=rgb(col,max=255),axes=TRUE, axisnames=F,ylim=c(0,2.8),cex.lab=1.7,cex.axis=1.7,beside=T,cex.names=1.7,space=0.5,ylab="Fold enrichment")
arrows(plot,t(enrich+enrichse),plot,t(enrich-enrichse), angle=90, code=3, length=0.03)
text(plot,par("usr")[3]-0.05,labels=names,srt=45,adj=1,xpd=TRUE,cex=1.7)
abline(h=1,lty=2,col=2)

dev.off()
########################################################
#            Compute the fold enrichment               #
######################################################## 
alter=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.anno.txt",head=T)
siganno=rowSums(alter[,2:dim(alter)[2]])/sum(rowSums(alter[,2:dim(alter)[2]]))
expanno=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/null_mprobes.txt",head=F)
enrich=siganno/apply(expanno,1,mean)
p_value=c()
for(i in 1:14){
    p=length(which(expanno[i,]>siganno[i]))/dim(expanno)[2]
    p_value=c(p_value,p)
}
output = data.frame(alter[,1],enrich,p_value)
output[order(output$enrich),]
enrichse = apply(siganno/expanno,1,sd)
#########################################################################################
#                                           End                                         #  
######################################################################################### 
