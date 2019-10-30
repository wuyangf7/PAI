# Partition the Hi-C loops from Jung et al. 2019 NG 
# to different confidence sets and
# select the PAIs overlap with either the true or 
# false sets based on the Hi-C p values
# and see the PAIs enrichment 
#############################################################
   # True and False Hi-C loops at different FDR threshold #
#############################################################
library(data.table); library(stringr)
setwd("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/")
PP = fread("GM.pp.txt",head=T,data.table=F,check.names=F)
PO = fread("GM.po.txt",head=T,data.table=F,check.names=F)
all = rbind(PP,PO)
####### bait
chr1 = str_split_fixed(all$frag1,":",2)[,1]
tmp1 = str_split_fixed(all$frag1,":",2)[,2]
start1 = str_split_fixed(tmp1,"-",2)[,1]
end1 = str_split_fixed(tmp1,"-",2)[,2]
####### target
chr2 = str_split_fixed(all$frag2,":",2)[,1]
tmp2 = str_split_fixed(all$frag2,":",2)[,2]
start2 = str_split_fixed(tmp2,"-",2)[,1]
end2 = str_split_fixed(tmp2,"-",2)[,2]
####### Partition
dat = data.frame(chr1,start1,end1,chr2,start2,end2,all$'-log10(result)')
thresh = c(0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.8,1,1.5,1.7,2:10)
for(t in 1:length(thresh)) {
	sigidx = which(dat[,7] >= thresh[t])
	true = dat[sigidx,]
	false = dat[-sigidx,]
	write.table(true,paste0("./Partition/true_",thresh[t],".bedpe"),row=F,col=F,sep="\t",qu=F)
	write.table(false,paste0("./Partition/false_",thresh[t],".bedpe"),row=F,col=F,sep="\t",qu=F)
}
############################################################# 
		 # sig. and non-sig. PAIs in bedpe format #
#############################################################
setwd("/shares/compbio/Group-Yang/y.wu/interaction/HiC/"); library(data.table); options(scipen=999)
all=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt",head=T,sep="\t",data.table=F)
all$Expo_Chr=paste0("chr",all$Expo_Chr); all$Outco_Chr=paste0("chr",all$Outco_Chr)
allbedpe=cbind(all$Expo_Chr,all$Expo_bp,all$Expo_bp,all$Outco_Chr,all$Outco_bp,all$Outco_bp)
sigidx = which(all$p_SMR <= 0.05/nrow(all) & all$p_HEIDI >= 0.01)
sig = allbedpe[sigidx,]; nonsig = allbedpe[-sigidx,]; 
sigidx = sample(1:nrow(allbedpe),nrow(all)/2)
#sigidx = sample(1:nrow(allbedpe),nrow(sig))
ransig = allbedpe[sigidx,]; rannonsig = allbedpe[-sigidx,] 
write.table(ransig,"/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/ransigDNAm_half.bedpe",row=F,,col=F,sep="\t",qu=F)
write.table(rannonsig,"/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/rannonsigDNAm_half.bedpe",row=F,,col=F,sep="\t",qu=F)
write.table(sig,"/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/sigDNAm.bedpe",row=F,,col=F,sep="\t",qu=F)
write.table(nonsig,"/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/nonsigDNAm.bedpe",row=F,,col=F,sep="\t",qu=F)
############################################################# 
				      # HiC overlap #
#############################################################
pgltools="/shares/compbio/Group-Yang/y.wu/interaction/HiC/pgltools-master/sh/pgltools"
sigDNAm="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/sigDNAm"
nonsigDNAm="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/nonsigDNAm"
ransigDNAm="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/ransigDNAm"
rannonsigDNAm="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/rannonsigDNAm"
junghicdir="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Partition/"
####### format Hi-C
thresh = c(0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.8,1,1.5,1.7,2:10)
for(t in 1:length(thresh)) {
	truecmd = paste0(pgltools," formatbedpe ",junghicdir,"true_",thresh[t],".bedpe > ",junghicdir,"true_",thresh[t],".pgl")
	falsecmd = paste0(pgltools," formatbedpe ",junghicdir,"false_",thresh[t],".bedpe > ",junghicdir,"false_",thresh[t],".pgl")
	system(truecmd)
	system(falsecmd)
}
####### overlap
thresh = c(0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.8,1,1.5,1.7,2:10)
for(t in 1:length(thresh)) {
	# tpcmd1 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",sigDNAm,".pgl -b ",junghicdir,"true_",thresh[t],".pgl -wa > ",junghicdir,"true_pos_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	# tncmd1 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",nonsigDNAm,".pgl -b ",junghicdir,"true_",thresh[t],".pgl -wa > ",junghicdir,"true_neg_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	# fpcmd1 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",sigDNAm,".pgl -b ",junghicdir,"false_",thresh[t],".pgl -wa > ",junghicdir,"false_pos_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	# fncmd1 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",nonsigDNAm,".pgl -b ",junghicdir,"false_",thresh[t],".pgl -wa > ",junghicdir,"false_neg_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	# system(tpcmd1); system(tncmd1);
	# system(fpcmd1); system(fncmd1);
	tpcmd2 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",ransigDNAm,".pgl -b ",junghicdir,"true_",thresh[t],".pgl -wa > ",junghicdir,"ran_true_pos_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	tncmd2 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",rannonsigDNAm,".pgl -b ",junghicdir,"true_",thresh[t],".pgl -wa > ",junghicdir,"ran_true_neg_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	fpcmd2 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",ransigDNAm,".pgl -b ",junghicdir,"false_",thresh[t],".pgl -wa > ",junghicdir,"ran_false_pos_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	fncmd2 = paste0("Rscript ~/bin/rsub.R '",pgltools," intersect -a ",rannonsigDNAm,".pgl -b ",junghicdir,"false_",thresh[t],".pgl -wa > ",junghicdir,"ran_false_neg_",thresh[t],"_overlap.pgl'"," @overlap %30 T1 W5")
	system(tpcmd2); system(tncmd2);
	system(fpcmd2); system(fncmd2); 
}
############################################################# 
		  # 	compute the fold enrichment 	 # 
#############################################################
junghicdir="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Partition/"
thresh = c(0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.8,1,1.5,1.7,2:10)
tpall1 = tnall1 = fpall1 = fnall1 = tpall2 = tnall2 = fpall2 = fnall2 = tpall3 = tnall3 = fpall3 = fnall3 = tpall4 = tnall4 = fpall4 = fnall4 = tpall5 = tnall5 = fpall5 = fnall5 = tpall6 = tnall6 = fpall6 = fnall6 = c()
for(t in 1:length(thresh)) {
	print(t)
	tp1 = system(paste0("cat ",junghicdir,"true_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tn1 = system(paste0("cat ",junghicdir,"true_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fp1 = system(paste0("cat ",junghicdir,"false_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fn1 = system(paste0("cat ",junghicdir,"false_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tp2 = system(paste0("cat ",junghicdir,"ran_true_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tn2 = system(paste0("cat ",junghicdir,"ran_true_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fp2 = system(paste0("cat ",junghicdir,"ran_false_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fn2 = system(paste0("cat ",junghicdir,"ran_false_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tp3 = system(paste0("cat ",junghicdir,"Corr_true_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tn3 = system(paste0("cat ",junghicdir,"Corr_true_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fp3 = system(paste0("cat ",junghicdir,"Corr_false_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fn3 = system(paste0("cat ",junghicdir,"Corr_false_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tp4 = system(paste0("cat ",junghicdir,"CACorr_true_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tn4 = system(paste0("cat ",junghicdir,"CACorr_true_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fp4 = system(paste0("cat ",junghicdir,"CACorr_false_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fn4 = system(paste0("cat ",junghicdir,"CACorr_false_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tp5 = system(paste0("cat ",junghicdir,"CASMR_true_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tn5 = system(paste0("cat ",junghicdir,"CASMR_true_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fp5 = system(paste0("cat ",junghicdir,"CASMR_false_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fn5 = system(paste0("cat ",junghicdir,"CASMR_false_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tp6 = system(paste0("cat ",junghicdir,"CAPHM_true_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tn6 = system(paste0("cat ",junghicdir,"CAPHM_true_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fp6 = system(paste0("cat ",junghicdir,"CAPHM_false_pos_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	fn6 = system(paste0("cat ",junghicdir,"CAPHM_false_neg_",thresh[t],"_overlap.pgl | wc -l"),intern=T)
	tpall1 = c(tpall1,tp1)
	tnall1 = c(tnall1,tn1)
	fpall1 = c(fpall1,fp1)
	fnall1 = c(fnall1,fn1)
	tpall2 = c(tpall2,tp2)
	tnall2 = c(tnall2,tn2)
	fpall2 = c(fpall2,fp2)
	fnall2 = c(fnall2,fn2)
	tpall3 = c(tpall3,tp3)
	tnall3 = c(tnall3,tn3)
	fpall3 = c(fpall3,fp3)
	fnall3 = c(fnall3,fn3)
	tpall4 = c(tpall4,tp4)
	tnall4 = c(tnall4,tn4)
	fpall4 = c(fpall4,fp4)
	fnall4 = c(fnall4,fn4)
	tpall5 = c(tpall5,tp5)
	tnall5 = c(tnall5,tn5)
	fpall5 = c(fpall5,fp5)
	fnall5 = c(fnall5,fn5)
	tpall6 = c(tpall6,tp6)
	tnall6 = c(tnall6,tn6)
	fpall6 = c(fpall6,fp6)
	fnall6 = c(fnall6,fn6)
}
tpall1=as.numeric(tpall1);tnall1=as.numeric(tnall1);fpall1=as.numeric(fpall1);fnall1=as.numeric(fnall1)
tpall2=as.numeric(tpall2);tnall2=as.numeric(tnall2);fpall2=as.numeric(fpall2);fnall2=as.numeric(fnall2)
tpall3=as.numeric(tpall3);tnall3=as.numeric(tnall3);fpall3=as.numeric(fpall3);fnall3=as.numeric(fnall3)
tpall4=as.numeric(tpall4);tnall4=as.numeric(tnall4);fpall4=as.numeric(fpall4);fnall4=as.numeric(fnall4)
tpall5=as.numeric(tpall5);tnall5=as.numeric(tnall5);fpall5=as.numeric(fpall5);fnall5=as.numeric(fnall5)
tpall6=as.numeric(tpall6);tnall6=as.numeric(tnall6);fpall6=as.numeric(fpall6);fnall6=as.numeric(fnall6)
orall1 = tpall1/fpall1/(tnall1/fnall1)
orall2 = tpall2/fpall2/(tnall2/fnall2)
orall3 = tpall3/fpall3/(tnall3/fnall3)
orall4 = tpall4/fpall4/(tnall4/fnall4)
orall5 = tpall5/fpall5/(tnall5/fnall5)
orall6 = tpall6/fpall6/(tnall6/fnall6)
ormat1 = matrix(0,4,length(thresh))
ormat2 = matrix(0,4,length(thresh))
ormat3 = matrix(0,4,length(thresh))
ormat4 = matrix(0,4,length(thresh))
ormat5 = matrix(0,4,length(thresh))
ormat6 = matrix(0,4,length(thresh))
ormat1[1,] = tpall1; ormat1[2,] = fpall1; ormat1[3,] = tnall1; ormat1[4,] = fnall1;
ormat2[1,] = tpall2; ormat2[2,] = fpall2; ormat2[3,] = tnall2; ormat2[4,] = fnall2;
ormat3[1,] = tpall3; ormat3[2,] = fpall3; ormat3[3,] = tnall3; ormat3[4,] = fnall3;
ormat4[1,] = tpall4; ormat4[2,] = fpall4; ormat4[3,] = tnall4; ormat4[4,] = fnall4;
ormat5[1,] = tpall5; ormat5[2,] = fpall5; ormat5[3,] = tnall5; ormat5[4,] = fnall5;
ormat6[1,] = tpall6; ormat6[2,] = fpall6; ormat6[3,] = tnall6; ormat6[4,] = fnall6;
############  confidence interval 
CIsiglow = CIsighigh = CIranlow = CIranhigh = CIcorrlow = CIcorrhigh = 
CICAcorrlow = CICAcorrhigh = CICASMRlow = CICASMRhigh = CICAPHMlow = CICAPHMhigh = c()
for(i in 1:length(thresh)) {
	lowtmp1 = fisher.test(matrix(ormat1[,i],nrow=2))$conf.int[1]
	hightmp1 = fisher.test(matrix(ormat1[,i],nrow=2))$conf.int[2]
	lowtmp2 = fisher.test(matrix(ormat2[,i],nrow=2))$conf.int[1]
	hightmp2 = fisher.test(matrix(ormat2[,i],nrow=2))$conf.int[2]
	lowtmp3 = fisher.test(matrix(ormat3[,i],nrow=2))$conf.int[1]
	hightmp3 = fisher.test(matrix(ormat3[,i],nrow=2))$conf.int[2]
	lowtmp4 = fisher.test(matrix(ormat4[,i],nrow=2))$conf.int[1]
	hightmp4 = fisher.test(matrix(ormat4[,i],nrow=2))$conf.int[2]
	lowtmp5 = fisher.test(matrix(ormat5[,i],nrow=2))$conf.int[1]
	hightmp5 = fisher.test(matrix(ormat5[,i],nrow=2))$conf.int[2]
	lowtmp6 = fisher.test(matrix(ormat6[,i],nrow=2))$conf.int[1]
	hightmp6 = fisher.test(matrix(ormat6[,i],nrow=2))$conf.int[2]
	CIsiglow = c(CIsiglow, lowtmp1)
	CIsighigh = c(CIsighigh, hightmp1)
	CIranlow = c(CIranlow, lowtmp2)
	CIranhigh = c(CIranhigh, hightmp2)
	CIcorrlow = c(CIcorrlow, lowtmp3)
	CIcorrhigh = c(CIcorrhigh, hightmp3)
	CICAcorrlow = c(CICAcorrlow, lowtmp4)
	CICAcorrhigh = c(CICAcorrhigh, hightmp4)
	CICASMRlow = c(CICASMRlow, lowtmp5)
	CICASMRhigh = c(CICASMRhigh, hightmp5)
	CICAPHMlow = c(CICAPHMlow, lowtmp6)
	CICAPHMhigh = c(CICAPHMhigh, hightmp6)
}
###########  plot 1 
pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/DNAm_overlap_JungHiC_threshold_enrichment.pdf")
plot(thresh[-c(1:5)],orall1[-c(1:5)],type="o",col="orange",pch=20,lty=1,lwd=2,xlab=expression(paste("-", log[10], "(PHi-C threshold ", italic(P), " value)", sep="")),ylab="Fold enrichment",main="",ylim=c(0,11),xlim=c(0,10),bty="n",cex.axis=1.5,cex.lab=1.5)
lines(thresh[-c(1:5)],orall2[-c(1:5)],type="o",col="blue",pch=18,lty=1,lwd=2)
legend("topleft",c(paste0("significant PAIs"),paste0("random PAIs")),lty=c(1,1),pch=c(20,18),col=c("orange","blue"),bty="n",cex=1.5)
abline(h=1,lty=2,col="red")
###########  plot 2 
pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/DNAm_overlap_JungHiC_threshold_enrichment_v5.pdf")
plot(thresh[-c(1:5,18:20)],orall1[-c(1:5,18:20)],type="o",col="orange",pch=19,lty=1,lwd=2,xlab=expression(paste("-", log[10], "(Hi-C threshold ", italic(P), " value)", sep="")),ylab="Fold enrichment",main="",ylim=c(0,13),xlim=c(0,8),bty="n",cex.axis=1.5,cex.lab=1.5)
segments(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)], lty=3,lwd=1.5,col="orange")
arrows(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)],lty=3, angle=90, code=3, length=0.05,lwd=2,col="orange")
points(thresh[-c(1:5,18:20)],orall2[-c(1:5,18:20)],type="o",col="blue",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="blue")
arrows(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="blue")
points(thresh[-c(1:5,18:20)],orall3[-c(1:5,18:20)],type="o",col="green",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="green")
arrows(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="green")
points(thresh[-c(1:5,18:20)],orall4[-c(1:5,18:20)],type="o",col="purple",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="purple")
arrows(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="purple")
points(thresh[-c(1:5,18:20)],orall5[-c(1:5,18:20)],type="o",col="mediumaquamarine",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="mediumaquamarine")
arrows(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="mediumaquamarine")
legend("topleft",c("Sig. PAIs using DNAm","Sig. PAIs using ATAC-seq","Correlation-based method using ATAC-seq","Correlation-based method using DNAm","random PAIs"),lty=c(1,1,1,1,1),pch=c(19,19,19,19,19),col=c("orange","mediumaquamarine","purple","green","blue"),bty="n",cex=1)
abline(h=1,lty=2,col="red") 
###########  plot 3 
pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/DNAm_overlap_JungHiC_threshold_enrichment_v6.pdf")
plot(thresh[-c(1:5,18:20)],orall1[-c(1:5,18:20)],type="o",col="orange",pch=19,lty=1,lwd=2,xlab=expression(paste("-", log[10], "(Hi-C threshold ", italic(P), " value)", sep="")),ylab="Fold enrichment",main="",ylim=c(0,13),xlim=c(0,8),bty="n",cex.axis=1.5,cex.lab=1.5)
segments(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)], lty=3,lwd=1.5,col="orange")
arrows(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)],lty=3, angle=90, code=3, length=0.05,lwd=2,col="orange")
points(thresh[-c(1:5,18:20)],orall2[-c(1:5,18:20)],type="o",col="blue",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="blue")
arrows(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="blue")
points(thresh[-c(1:5,18:20)],orall3[-c(1:5,18:20)],type="o",col="green",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="green")
arrows(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="green")
points(thresh[-c(1:5,18:20)],orall4[-c(1:5,18:20)],type="o",col="purple",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="purple")
arrows(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="purple")
points(thresh[-c(1:5,18:20)],orall5[-c(1:5,18:20)],type="o",col="mediumaquamarine",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="mediumaquamarine")
arrows(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="mediumaquamarine")
points(thresh[-c(1:5,18:20)],orall6[-c(1:5,18:20)],type="o",col="black",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICAPHMhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAPHMlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="black")
arrows(thresh[-c(1:5,18:20)],CICAPHMhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAPHMlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="black")
legend("topleft",c("Sig. PAIs using DNAm","Sig. PAIs using ATAC-seq","PHM using ATAC-seq","Correlation-based method using ATAC-seq","Correlation-based method using DNAm","random PAIs"),lty=c(1,1,1,1,1,1),pch=c(19,19,19,19,19,19),col=c("orange","mediumaquamarine","black","purple","green","blue"),bty="n",cex=1)
abline(h=1,lty=2,col="red") 
#############################################################

#############################################################
	 #combined plot wiht CHIA-PET enrich (08.10.2019)#
	 			# include PHM results #
#############################################################
CHIAPETobsv=2315;
CHIAPETnull=read.table("/shares/compbio/Group-Yang/y.wu/interaction/HiC/enrich/res/nullDNAm_001_distM_rmMHC_CHIAPET.txt",head=F)
CHIAPET=c(CHIAPETnull[,1],CHIAPETobsv)
pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/DNAm_overlap_JungHiC_threshold_enrichment_plus_ChIA-PET.pdf",height=8,width=16)
par(mar=c(5,5,3,2),mfrow=c(1,2))
h=hist(CHIAPETnull[,1],breaks=50,col=c(rep("darkolivegreen1",215)),xlim=c(1400,2500),xlab="Number of overlaps with loops",main="ChIA-PET (POLR2A)",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
abline(v=CHIAPETobsv,lty=2,lwd=2,col="red"); 
mtext("a",side=3,line=0.5,adj=-0.15,cex=3.5,col="black",font=1)
plot(thresh[-c(1:5,18:20)],orall1[-c(1:5,18:20)],type="o",col="orange",pch=19,lty=1,lwd=2,xlab=expression(paste("-", log[10], "(PCHi-C ", italic(P), " value threshold)", sep="")),ylab="Fold enrichment",main="",ylim=c(0,13),xlim=c(0,8),bty="n",cex.axis=1.5,cex.lab=1.5)
segments(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)], lty=3,lwd=1.5,col="orange")
arrows(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)],lty=3, angle=90, code=3, length=0.05,lwd=2,col="orange")
points(thresh[-c(1:5,18:20)],orall2[-c(1:5,18:20)],type="o",col="blue",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="blue")
arrows(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="blue")
points(thresh[-c(1:5,18:20)],orall3[-c(1:5,18:20)],type="o",col="green",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="green")
arrows(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="green")
points(thresh[-c(1:5,18:20)],orall4[-c(1:5,18:20)],type="o",col="purple",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="purple")
arrows(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="purple")
points(thresh[-c(1:5,18:20)],orall5[-c(1:5,18:20)],type="o",col="mediumaquamarine",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="mediumaquamarine")
arrows(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="mediumaquamarine")
points(thresh[-c(1:5,18:20)],orall6[-c(1:5,18:20)],type="o",col="black",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICAPHMhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAPHMlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="black")
arrows(thresh[-c(1:5,18:20)],CICAPHMhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAPHMlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="black")
legend("topleft",c("Sig. PAIs using DNAm","Sig. PAIs using ATAC-seq","PHM using ATAC-seq","Correlation-based method using ATAC-seq","Correlation-based method using DNAm","random PAIs"),lty=c(1,1,1,1,1,1),pch=c(19,19,19,19,19,19),col=c("orange","mediumaquamarine","black","purple","green","blue"),bty="n",cex=1.5)
abline(h=1,lty=2,col="red")
mtext("b",side=3,line=0.1,adj=-0.15,cex=3.5,col="black",font=1) 
dev.off()
#############################################################
	 #combined plot wiht CHIA-PET enrich (09.10.2019)# 
	 			  # exclude PHM results #
#############################################################
CHIAPETobsv=2315;
CHIAPETnull=read.table("/shares/compbio/Group-Yang/y.wu/interaction/HiC/enrich/res/nullDNAm_001_distM_rmMHC_CHIAPET.txt",head=F)
CHIAPET=c(CHIAPETnull[,1],CHIAPETobsv)
pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/DNAm_overlap_JungHiC_threshold_enrichment_plus_ChIA-PET_v1.pdf",height=8,width=16)
par(mar=c(5,5,3,2),mfrow=c(1,2))
h=hist(CHIAPETnull[,1],breaks=50,col=c(rep("darkolivegreen1",215)),xlim=c(1400,2500),xlab="Number of overlaps with loops",main="ChIA-PET (POLR2A)",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
abline(v=CHIAPETobsv,lty=2,lwd=2,col="red"); 
mtext("a",side=3,line=0.5,adj=-0.15,cex=3.5,col="black",font=1)
plot(thresh[-c(1:5,18:20)],orall1[-c(1:5,18:20)],type="o",col="orange",pch=19,lty=1,lwd=2,xlab=expression(paste("-", log[10], "(Hi-C threshold ", italic(P), " value)", sep="")),ylab="Fold enrichment",main="",ylim=c(0,13),xlim=c(0,8),bty="n",cex.axis=1.5,cex.lab=1.5)
segments(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)], lty=3,lwd=1.5,col="orange")
arrows(thresh[-c(1:5,18:20)],CIsighigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIsiglow[-c(1:5,18:20)],lty=3, angle=90, code=3, length=0.05,lwd=2,col="orange")
points(thresh[-c(1:5,18:20)],orall2[-c(1:5,18:20)],type="o",col="blue",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="blue")
arrows(thresh[-c(1:5,18:20)],CIranhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIranlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="blue")
points(thresh[-c(1:5,18:20)],orall3[-c(1:5,18:20)],type="o",col="green",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="green")
arrows(thresh[-c(1:5,18:20)],CIcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CIcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="green")
points(thresh[-c(1:5,18:20)],orall4[-c(1:5,18:20)],type="o",col="purple",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="purple")
arrows(thresh[-c(1:5,18:20)],CICAcorrhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICAcorrlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="purple")
points(thresh[-c(1:5,18:20)],orall5[-c(1:5,18:20)],type="o",col="mediumaquamarine",pch=19,lty=1,lwd=2)
segments(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], lty=3,lwd=1.5,col="mediumaquamarine")
arrows(thresh[-c(1:5,18:20)],CICASMRhigh[-c(1:5,18:20)],thresh[-c(1:5,18:20)],CICASMRlow[-c(1:5,18:20)], angle=90, code=3, length=0.05,lty=3,lwd=2,col="mediumaquamarine")
legend("topleft",c("Sig. PAIs using DNAm","Sig. PAIs using ATAC-seq","Correlation-based method using ATAC-seq","Correlation-based method using DNAm","random PAIs"),lty=c(1,1,1,1,1),pch=c(19,19,19,19,19),col=c("orange","mediumaquamarine","purple","green","blue"),bty="n",cex=1)
abline(h=1,lty=2,col="red")
mtext("b",side=3,line=0.1,adj=-0.15,cex=3.5,col="black",font=1) 
dev.off()

#############################################################
	 #combined plot wiht CHIA-PET enrich (23.10.2019)#
	 		# include PHM results and bar plot #
#############################################################
CHIAPETobsv=2315;HiCobsv=130;
CHIAPETnull=read.table("/shares/compbio/Group-Yang/y.wu/interaction/HiC/enrich/res/nullDNAm_001_distM_rmMHC_CHIAPET.txt",head=F)
HiCnull=read.table("/shares/compbio/Group-Yang/y.wu/interaction/HiC/enrich/res/nullDNAm_001_distM_rmMHCdup_HiC.txt",head=F)
CHIAPET=c(CHIAPETnull[,1],CHIAPETobsv)
HiC=c(HiCnull[,1],HiCobsv)

pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/DNAm_overlap_JungHiC_threshold_enrichment_plus_ChIA-PET_plus_HiC.pdf",height=9,width=9)
par(mar=c(5,5,3,2))
layout(matrix(c(1,2,3,3),nrow=2,byrow=T))
#####		Hi-C loops from Rao et al.
h=hist(HiCnull[,1],breaks=50,col=c(rep("peru",215)),xlim=c(60,140),xlab="Number of overlaps with loops",main="Hi-C (Rao et al.)",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
abline(v=HiCobsv,lty=2,lwd=2,col="red"); #mtext("b",side=3,line=0.5,adj=-0.15,cex=1.8,col="black",font=1)
mtext("a",side=3,line=0.5,adj=-0.2,cex=2.5,col="black",font=1)
#####		ChIA-PET loops
h=hist(CHIAPETnull[,1],breaks=50,col=c(rep("darkolivegreen1",215)),xlim=c(1400,2500),xlab="Number of overlaps with loops",main="ChIA-PET (POLR2A)",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
abline(v=CHIAPETobsv,lty=2,lwd=2,col="red"); 
mtext("b",side=3,line=0.5,adj=-0.2,cex=2.5,col="black",font=1)
#####		PCHi-C loops from Jung et al.
data=cbind(orall1,orall5,orall6,orall4,orall3,orall2); 
dataCIlow=cbind(CIsiglow,CICASMRlow,CICAPHMlow,CICAcorrlow,CIcorrlow,CIranlow);
dataCIhigh=cbind(CIsighigh,CICASMRhigh,CICAPHMhigh,CICAcorrhigh,CIcorrhigh,CIranhigh);
thresh = c(0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.8,1,1.5,1.7,2:10)
#####
rowidx = c(9,12,13,14,15,16,17); data = data[rowidx,]; dataCIlow = dataCIlow[rowidx,]
dataCIhigh = dataCIhigh[rowidx,]; thresh = thresh[rowidx]
space = c(c(rep(0,6),1.5),rep(c(rep(0,5),1.5),length(rowidx)-2),rep(0,5))
#####
plot<-barplot(t(data),col=c("orange","mediumaquamarine","black","purple","green","blue"),space=space,axes=TRUE, axisnames=F,ylim=c(0,13),cex.lab=1.5,beside=T,cex.names=1.5,ylab="Fold enrichment",xlab=expression(paste("-", log[10], "(PCHi-C ", italic(P), " value threshold)", sep="")),cex.axis=1.5,main="",cex.main=1.5)
legend("topleft",c("SMR & HEIDI using DNAm","SMR & HEIDI using ATAC-seq","PHM using ATAC-seq","Correlation-based method using ATAC-seq","Correlation-based method using DNAm","Random PAIs"),cex=1.5,fill=c("orange","mediumaquamarine","black","purple","green","blue"),bty="n")
names = thresh
axis(1, at=c(3, 10.5, 18, 25.5, 33, 40.5,48), labels=names, srt=45,las=1, cex.axis=1.5)
segments(plot,t(dataCIhigh),plot,t(dataCIlow), lwd=1.5)
arrows(plot,t(dataCIhigh),plot,t(dataCIlow), angle=90, code=3, length=0.03)
abline(h=1,lty=2,col="red")
mtext("c",side=3,line=0.5,adj=-0.1,cex=2.5,col="black",font=1)
dev.off()
#############################################################
#############################################################