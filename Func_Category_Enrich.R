############################################################################
#       Compute the enrichment of PIDSs in 14 functional categories        # 
#                       from the Roadmap project                           # 
#                   1. Estimate the fold enrichment                        # 
#       2. caculate se based on 100 null enrichment distribution           # 
#                        3. summarize and plot                              # 
############################################################################
#########################################################################################
#                         1. estimate the fold enrichment                               #  
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
for(i in 1:500){
    print(i)
    smplmprobes=c()
    for(j in 1:nvarbin)
    {
        set.seed(i)
        smplmprobes=c(smplmprobes,sample(splsets[[j]],nperitvl[j],replace = FALSE))
    }
    outfile=paste0("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/sampledProbes/null_mprobes_",i,".bim")
    write.table(allDNAm[smplmprobes,c(1:3)],outfile,row=F,col=F,sep="\t",quote=F)
} 
######          observed PIDSs formatting
write.table(sigDNAm,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.bim",row=F,col=F,sep="\t",quote=F)
########################################################
#            Annotation for null samples               # 
########################################################
anno="/shares/compbio/Group-Yang/y.wu/omics/enrichment/script/anno14imb"
for((i=1; i<=500; i++))
do
    bim="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/sampledProbes/null_mprobes_"$i".bim"
    out="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/jobs/null_mprobes_"$i
    log="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/log/null_mprobes_"$i".log"
    rsub "$anno --cat 1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
    --bim $bim --out $out >$log 2>&1" @anno %3 T1 W1
done
########################################################
#                Annotation for PIDSs                  # 
########################################################
anno="/shares/compbio/Group-Yang/y.wu/omics/enrichment/script/anno14imb"
$anno --cat 1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
--bim /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.bim \
--out /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr \
> /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.log 2>&1
########################################################
#            Annotation for all DNAm probes            # 
########################################################
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
simu=500
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
#                   2. estimated se from 100 null enrichment analysis                   #  
#########################################################################################
sigfile = "/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt"
allfile = "/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt"
sigDNAm = unique(read.table(sigfile,header=T,sep="\t")[,c(6,5,8)])
allDNAm = unique(fread(allfile,head=T,data.table=F)[,c(6,5,8)])
mevar = read.table("/shares/compbio/Group-Yang/y.wu/interaction/data/lbc.vp.txt",head=F)
indx = match(allDNAm$Outco_ID,mevar$V1,nomatch=0)
allDNAm$var = mevar$V2[indx]
options(scipen=999)
sigindx = match(sigDNAm$Outco_ID,allDNAm$Outco_ID,nomatch=0)
sigvar = allDNAm$var[sigindx]
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
for(rep in 1:101){
    outdir=paste("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/sampledProbes/rep",rep,"/",sep="")
    dir.create(outdir,showWarnings = FALSE)
    for(i in 1:100){
        print(i)
        smplmprobes=c()
        for(j in 1:nvarbin)
        {
            smplmprobes=c(smplmprobes,sample(splsets[[j]],nperitvl[j],replace = FALSE))
        }
        outfile=paste0(outdir,"null_mprobes_",i,".bim")
        write.table(allDNAm[smplmprobes,c(1:3)],outfile,row=F,col=F,sep="\t",quote=F)
    }
}
#####################################################  
#                   Annotation                      #
#####################################################
anno="/shares/compbio/Group-Yang/y.wu/omics/enrichment/script/anno14imb"
for((r=1; r<=101; r++))
do
for((i=1; i<=100; i++))
do
bim="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/sampledProbes/rep"$r"/null_mprobes_"$i".bim"
out="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/output/rep"$r"/null_mprobes_"$i
log="/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/logfile/rep"$r"/null_mprobes_"$i".log"
rsub "$anno --cat 1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
--bim $bim --out $out >$log 2>&1" @enrich %3 T1 W3
done
done
#####################################################
#           Summarize enrich each replicate         #
#####################################################
for((r=1; r<=101; r++))
do
rsub "Rscript /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/summarize.r $r \
> /shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/logfile/${r}.log 2>&1" @sum %5 T1 W1
done
#####################################################
#                    Estimate se                    #
#####################################################
library(data.table)
nullrep=matrix(0,14,100)
for(r in 1:100){
print(r)
alter=read.table(paste0("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/output/rep1/null_mprobes_",r,".anno.txt"),head=T)
nullsig=rowSums(alter[,2:dim(alter)[2]])/sum(rowSums(alter[,2:dim(alter)[2]]))
nullexp=read.table(paste0("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/se/output/rep",r+1,"/null_mprobes_all_anno.txt"),head=F)
nullrep[,r]=nullsig/apply(nullexp,1,mean)
}
sigse=apply(nullrep,1,sd)

write.table(nullrep,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/mprobes_null_fold.txt",row=F,col=F,sep="\t",qu=F)
write.table(sigse,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/mprobes_enrich_se.txt",row=F,col=F,sep="\t",qu=F)
#########################################################################################
#                                 3. Summarize and plot                                 #  
######################################################################################### 
alter=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.anno.txt",head=T)
siganno=rowSums(alter[,2:dim(alter)[2]])/sum(rowSums(alter[,2:dim(alter)[2]]))
expanno=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/null_mprobes.txt",head=F)
enrich=siganno/apply(expanno,1,mean)
enrichse=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/mprobes_enrich_se.txt",head=F)
tvalue=as.matrix(abs(enrich-1)/enrichse)
pvalue=pt(tvalue,dim(alter)[2]-1,lower.tail=F)
data.frame(alter[,1],enrich,pvalue)

library(stringr)
colfile=unique(read.delim("/gpfs/gpfs01/polaris/Q0286/uqywu16/chromstate/annotation_25_imputed12marks.txt",head=T,sep="\t")$COLOR.CODE)
col=str_split_fixed(colfile,",",3)

pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/new/sigmprobes_enrich_smr.pdf",height=8,width=15)
par(mar=c(5,5,3,2),mfrow=c(1,2))
library(data.table)
obsvfile=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.anno.txt",head=T)
obsv=rowSums(obsvfile[,2:dim(obsvfile)[2]])/sum(rowSums(obsvfile[,2:dim(obsvfile)[2]]))
expfile=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/allprobes.anno.txt",head=T,data.table=F)
exp=rowSums(expfile[,2:dim(expfile)[2]])/sum(rowSums(expfile[,2:dim(expfile)[2]]))
data=cbind(obsv,exp)
names=c("TssA","Prom","Tx","TxWk","TxEn","EnhA","EnhW","DNase","ZNF/Rpts","Het","PromP","PromBiv","ReprPC","Quies")
plot<-barplot(t(data),col=c("cornflowerblue","forestgreen"),axes=TRUE, axisnames=F,ylim=c(0,0.35),cex.lab=1.7,cex.axis=1.7,beside=T,cex.names=1.7,ylab="Proportion of probes") 
legend("topright",c("Sig. DNAm","All DNAm"),cex=1.7,fill=c("cornflowerblue","forestgreen"),bty="n")
text(apply(plot,2,mean),par("usr")[3]-0.005,labels=names,srt=45,adj=1,xpd=TRUE,cex=1.7)
names=c("TssA","Prom","Tx","TxWk","TxEn","EnhA","EnhW","DNase","ZNF/Rpts","Het","PromP","PromBiv","ReprPC","Quies")   
plot<-barplot(enrich,col=rgb(col,max=255),axes=TRUE, axisnames=F,ylim=c(0,2.8),cex.lab=1.7,cex.axis=1.7,beside=T,cex.names=1.7,space=0.5,ylab="Fold enrichment")
arrows(plot,t(enrich+enrichse),plot,t(enrich-enrichse), angle=90, code=3, length=0.03)
text(plot,par("usr")[3]-0.05,labels=names,srt=45,adj=1,xpd=TRUE,cex=1.7)
abline(h=1,lty=2,col=2)
dev.off()
####################################################################
                # compute the enrichment P value # 
#################################################################### 
alter=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/sig_mprobes_smr.anno.txt",head=T)
expanno=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/null_mprobes.txt",head=F)
nullrep=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/FuncEnrich/res/mprobes_null_fold.txt",head=F)
siganno=rowSums(alter[,2:dim(alter)[2]])/sum(rowSums(alter[,2:dim(alter)[2]]))
enrich=siganno/apply(expanno,1,mean)

p_value=c()
for(i in 1:14){
    if(enrich[i] > 1) {p=length(which(nullrep[i,]>enrich[i]))/dim(nullrep)[2]
        } else {p=length(which(nullrep[i,]<enrich[i]))/dim(nullrep)[2]}
    p_value=c(p_value,p)
}
output = data.frame(alter[,1],enrich,p_value)
#########################################################################################
#                                           End                                         #  
######################################################################################### 