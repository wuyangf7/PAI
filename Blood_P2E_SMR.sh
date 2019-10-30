########################################## 
			#DNAm interactions
##########################################
COHORT="Blood"
SMR="/shares/compbio/Group-Yang/uqfzhan7/bin/smr_linux"
#SMR="/shares/compbio/Group-Yang/y.wu/bin/smr_linux"
METHY="/shares/compbio/PCTG/meQTL/LBC_BSGS/meta_meqtl_chr"
CAGE="/shares/compbio/Group-Yang/eQTL/CAGE/CAGE.cis.chr"
ProProbes="/shares/compbio/Group-Yang/y.wu/allSMR/m2msmr/P2E/promoter_mprobes.list"
HEIDIEXPO="/shares/compbio/Group-Yang/y.wu/interaction/P2E/res/DNAm_expo_pass_smr.list"
HEIDIOUTCO="/shares/compbio/Group-Yang/y.wu/interaction/P2E/res/DNAm_outco_pass_smr.list"
REFERENCE="/shares/compbio/Group-Yang/reference_data/hrs/hrs_1kg_hwe1e-6"
INDI="/shares/compbio/Group-Yang/uqzzhu1/dataset/SepIDs/hrs_hm3_m01_hwe_1e-6_u05_geno.grm.id"
HRSSNPLIST="/shares/compbio/Group-Yang/t.qi/data/SepIDs/hrs_1kg_impRsq_0.3.snplist"
TARGET="/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/add/target.list"
EXPROBES="/shares/compbio/Group-Yang/y.wu/omic2/QC/excl_mhc_tqtl_mprbs.list"
EXSNPS="/shares/compbio/Group-Yang/y.wu/SMRdata/MHC_snps.list"
GENELIST="/shares/compbio/Group-Yang/y.wu/database/refseq_isoform_SMR.txt"
###########################################  
			 #P2E (including)  
###########################################
OUTPUT="/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/"
for((i=1;i<=22;i++))
do
rsub "${SMR} --beqtl-summary ${METHY}${i} --extract-exposure-probe ${ProProbes} \
 --beqtl-summary ${METHY}${i} \
 --bfile ${REFERENCE}_chr${i} --exclude-snp ${EXSNPS} --extract-snp ${HRSSNPLIST} --keep ${INDI} \
 --maf 0.01 --cis-wind 2000 --out ${OUTPUT}${COHORT}_m2msmr_chr${i} --thread-num 20 \
 > ${OUTPUT}${COHORT}_m2msmr_chr${i}.log 2>&1" @interaction %4 T7 W150
done
###########################################  
			 #P2E (additional)  
###########################################
OUTPUT="/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/add/"
i=2
rsub "${SMR} --beqtl-summary ${METHY}${i} --extract-exposure-probe ${TARGET} \
 --beqtl-summary ${METHY}${i} \
 --bfile ${REFERENCE}_chr${i} --exclude-snp ${EXSNPS} --extract-snp ${HRSSNPLIST} --keep ${INDI} \
 --maf 0.01 --cis-wind 2000 --out ${OUTPUT}${COHORT}_m2msmr_chr${i} --thread-num 20 \
 > ${OUTPUT}${COHORT}_m2msmr_chr${i}.log 2>&1" @interaction %4 T7 W150
############################################ 
			#exclude same promoter 
############################################
promp=read.table("/shares/compbio/Group-Yang/y.wu/interaction/data/methylation_promoters.txt",head=F)
megeo=read.table("/shares/compbio/Group-Yang/y.wu/database/GPL16304-47833.txt",skip=22,head=T,sep="\t")
for (i in 1:22){
print(i)
dat=read.delim(paste("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/Blood_m2msmr_chr",i,".smr",sep=""),head=T,sep="\t")
indx=match(dat$Expo_ID,promp$V2,nomatch=0)
tmp=which(dat$Outco_bp>promp$V5[indx] & dat$Outco_bp<=promp$V6[indx])
data=dat[-tmp,]
expoindx=match(data$Expo_ID,megeo$ID,nomatch=0)
outcindx=match(data$Outco_ID,megeo$ID,nomatch=0)
data$Expo_Gene=megeo$Closest_TSS_gene_name[expoindx]
data$Outco_Gene=megeo$Closest_TSS_gene_name[outcindx]
write.table(data,paste("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_chr",i,".smr",sep=""),row=F,qu=F,sep="\t")}
############################################# 
				#summarize  
#############################################
setwd("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/")
datall=c()
for(i in 1:22) {
	print(i)
	file=paste0("Blood_m2msmr_chr",i,".smr")
	if(file.exists(file)==1){
	data=read.table(file,head=T,sep="\t",colClass=c(rep(c("character","numeric"),5),"numeric",rep("character",2),rep("numeric",12)))
	datall=rbind(datall,data)
	}
}
smrindx=which(datall$p_SMR<=(0.05/dim(datall)[1]))
heidiindx=which(datall[smrindx,]$p_HEIDI>=0.01)
write.table(datall,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all.txt",row=F,sep="\t",qu=F)
write.table(datall[smrindx,],"/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_pass_smr.txt",row=F,sep="\t",qu=F)
write.table(datall[smrindx[heidiindx],],"/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_pass_heidi.txt",row=F,sep="\t",qu=F)
############################################# 
add=read.table("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/add/Blood_m2msmr_chr2.smr",head=T)
all=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi.txt",head=T,sep="\t")
sigadd = which(add$p_SMR <= 1.57e-9& add$p_HEIDI >= 0.01)
all = rbind(all,add[sigadd,])
write.table(all,"/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi.txt",row=F,sep="\t",qu=F)
############################################# 
smr=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t",stringsAsFactors=F)
megeo=read.table("/shares/compbio/Group-Yang/y.wu/database/GPL16304-47833.txt",skip=22,head=T,sep="\t")
expoindx=match(smr$Expo_ID,megeo$ID,nomatch=0)
outcindx=match(smr$Outco_ID,megeo$ID,nomatch=0)
smr$Expo_Gene=megeo$Closest_TSS_gene_name[expoindx]
smr$Outco_Gene=megeo$Closest_TSS_gene_name[outcindx]
write.table(smr,"/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",row=F,qu=F,sep="\t")
############################################# 
				#remove MHC  
#############################################
library(data.table)
all = fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all.txt",head=T,data.table=F)
smr = read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi.txt",head=T,sep="\t",stringsAsFactors=F)
mhcStart = 26056121; mhcEnd = 33699377
smrmhc = which( (smr$Expo_Chr==6 & smr$Expo_bp>=mhcStart & smr$Expo_bp<=mhcEnd) | (smr$Outco_Chr==6 & smr$Outco_bp>=mhcStart & smr$Outco_bp<=mhcEnd) )
allmhc = which( (all$Expo_Chr==6 & all$Expo_bp>=mhcStart & all$Expo_bp<=mhcEnd) | (all$Outco_Chr==6 & all$Outco_bp>=mhcStart & all$Outco_bp<=mhcEnd) )
smr = smr[-smrmhc,]
all = all[-allmhc,]
outfile1 = "/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt"
outfile2 = "/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc.txt"
write.table(smr,outfile1,row=F,sep="\t",qu=F)
write.table(all,outfile2,row=F,sep="\t",qu=F)
############################################# 
				#Unique pairs  
#############################################
library(data.table)
all = fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc.txt",head=T,data.table=F)
smr = read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t",stringsAsFactors=F)
allpair1 = paste(all$Expo_ID, all$Outco_ID, sep="-")
smrpair1 = paste(smr$Expo_ID, smr$Outco_ID, sep="-")
allpair2 = paste(all$Outco_ID, all$Expo_ID, sep="-")
smrpair2 = paste(smr$Outco_ID, smr$Expo_ID, sep="-")
smrdupidx = match(smrpair2,smrpair1,nomatch=0); 
alldupidx = match(allpair2,allpair1,nomatch=0); 

smrrmidx = c(); index1=smrdupidx[which(smrdupidx!=0)]; index2=which(smrdupidx!=0)
for(i in 1:length(index1)) {
	idx1 = index1[i]; 	idx2 = index2[i]
	if((smr[idx1,]$b_Expo/smr[idx1,]$se_Expo)^2 >= (smr[idx2,]$b_Expo/smr[idx2,]$se_Expo)^2) {
		smrrmidx = c(smrrmidx,idx2)
	}
}
allrmidx = c(); index1=alldupidx[which(alldupidx!=0)]; index2=which(alldupidx!=0)
for(i in 1:length(index1)) {
	idx1 = index1[i]; 	idx2 = index2[i]
	if((all[idx1,]$b_Expo/all[idx1,]$se_Expo)^2 >= (all[idx2,]$b_Expo/all[idx2,]$se_Expo)^2) {
		allrmidx = c(allrmidx,idx2)
	}
}
smr = smr[-smrrmidx,]; all = all[-allrmidx,]
outfile1 = "/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt"
outfile2 = "/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt"
write.table(smr,outfile1,row=F,sep="\t",qu=F)
write.table(all,outfile2,row=F,sep="\t",qu=F)
############################################# 
			#Methy distribution#  
#############################################
library(data.table)
all = fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt",head=T,data.table=F)
smr = read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt",head=T,sep="\t")
load("/shares/compbio/Group-Yang/uqfzhan7/Data/LBC/Methylation/Beta_1342_wave1.RObject")
avermethy = apply(dat,1,mean)

sigDNAm = unique(c(as.character(smr$Expo_ID))); allDNAm = unique(c(as.character(all$Expo_ID)))
sigidx = match(sigDNAm, names(avermethy),nomatch=0); allidx = match(allDNAm, names(avermethy),nomatch=0)
pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/new/DNAm_distr_submit.pdf",height=5,width=10)
par(mar=c(5,5,3,2),mfrow=c(1,2))
hist(avermethy,breaks=80,col="blue",xlab="Mean methylation beta values",main="All DNAm probes",cex.axis=1,cex.lab=1,cex.main=1)
mtext("a",side=3,line=0.5,adj=-0.15,cex=1.8,col="black",font=1)
hist(avermethy[allidx],breaks=80,col="orange",xlab="Mean methylation beta values",main="Exposure DNAm probes",cex.axis=1,cex.lab=1,cex.main=1)
mtext("b",side=3,line=0.5,adj=-0.15,cex=1.8,col="black",font=1)
############################################# 
				#Shiny APP  
#############################################
library(shiny)
shiny::runApp("~/Desktop/DNAm/M2Mdb")
############################################# 
			   #Link with GWAS  
#############################################
OUTPUT="/shares/compbio/Group-Yang/y.wu/omics/m2tsmr/t2dnew/"
t2d="/shares/compbio/Group-Yang/a.xue/data/meta_DGU_whole_sum_for_SMR.txt"
ea="/shares/compbio/Group-Yang/reference_data/brain_trait_gwas/GWAS_EA3_excl23andMe.raw"
ibd="/shares/compbio/Group-Yang/y.wu/SMRdata/public_meta/ibd_iibdgc_2015.raw"
for((j=1;j<=22;j++))
do
rsub "${SMR}  --bfile ${REFERENCE}_chr${j} --beqtl-summary ${METHY}${j}\
 --gwas-summary $t2d --exclude-probe ${EXPROBES} --extract-snp ${HRSSNPLIST} --keep ${INDI}\
 --maf 0.01 --out ${OUTPUT}/t2dnew_m2tsmr_chr${j}\
 --thread-num 20 > ${OUTPUT}/t2dnew_m2tsmr_chr${j}.log 2>&1" @m2tsmr %10 T1 W20
done
OUTPUT="/shares/compbio/Group-Yang/y.wu/omics/g2tsmr/ibd"
for((j=1;j<=22;j++))
do
rsub "${SMR}  --bfile ${REFERENCE}_chr${j} --beqtl-summary ${CAGE}${j}\
 --gwas-summary $ibd --extract-snp ${HRSSNPLIST} --keep ${INDI}\
 --maf 0.01 --out ${OUTPUT}/ibd_g2tsmr_chr${j}\
 --thread-num 1 > ${OUTPUT}/ibd_g2tsmr_chr${j}.log 2>&1" @g2tsmr %10 T1 W20
done
############################################# 
			   #Match with GWAS  
#############################################
trait=c("height","bmi","whradjbmi","hdl","ldl","tg","ea","ra","scz","cad","t2d","cd","uc","ad","ibd")
m2t_heidi="/shares/compbio/Group-Yang/y.wu/omics/m2tsmr/pass_heidi/"

m2t=c()
for(t in 1:length(trait)){
m2t_tmp=read.table(paste(m2t_heidi,trait[t],"_m2tsmr_pass_heidi.smr",sep=""),head=T)
m2t_tmp=cbind(rep(trait[t],dim(m2t_tmp)[1]),m2t_tmp)
m2t=rbind(m2t,m2t_tmp)
}
colnames(m2t)[1]="Trait"
mhcStart = 26056121; mhcEnd = 33699377
m2tmhc = which( (m2t$ProbeChr==6 & m2t$Probe_bp>=mhcStart & m2t$Probe_bp<=mhcEnd) )
m2t = m2t[-m2tmhc,]
#############################################
#m2t=read.table("/shares/compbio/Group-Yang/y.wu/omic2/result/m2tsmr_pass_heidi.txt",head=T)
m2m=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc.txt",head=T,sep="\t")
megeo=read.table("/shares/compbio/Group-Yang/y.wu/database/GPL16304-47833.txt",skip=22,head=T,sep="\t")
DNAm=unique(c(as.character(m2m$Expo_ID),as.character(m2m$Outco_ID)))
idx=match(as.character(m2t$probeID),DNAm,nomatch=0)
m2tPAI=m2t[which(idx!=0),]
m2tPAI$Gene=megeo$Closest_TSS_gene_name[match(m2tPAI$probeID,megeo$ID)]
write.table(m2tPAI,"/shares/compbio/Group-Yang/y.wu/interaction/res/Trait_blood_DNAm_heidi.txt",row=F,qu=F,sep="\t")
#############################################
m2t=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Trait_blood_DNAm_heidi.txt",head=T,sep="\t")
###      g2tsmr      ###
trait=c("height","bmi","whradjbmi","hdl","ldl","tg","ea","ra","scz","cad","t2d","cd","uc","ad","ibd")
epi = read.table("/shares/compbio/Group-Yang/eQTL/CAGE/CAGE.cis.probes.epi",head=F)
g2t_heidi="/shares/compbio/Group-Yang/y.wu/omics/g2tsmr/pass_heidi/"

g2t=c()
for(t in 1:length(trait)){
	g2t_tmp=read.table(paste(g2t_heidi,trait[t],"_g2tsmr_cage_pass_heidi.smr",sep=""),head=T)
	g2t_tmp=cbind(rep(trait[t],dim(g2t_tmp)[1]),g2t_tmp)
	g2t=rbind(g2t,g2t_tmp)
}
colnames(g2t)[1]="Trait"
pairg2t = unique(paste(g2t$Trait,g2t$Gene,sep="-"))
pairm2t = unique(paste(m2t$Trait,m2t$Gene,sep="-"))
idx = match(pairm2t, pairg2t, nomatch=0)
length(which(idx!=0))/length(idx)
unmatch = unique(str_split_fixed(pairm2t[which(idx==0)],"-",2)[,2])
gidx = match(unmatch, epi$V5, nomatch=0)
length(which(gidx==0))/length(gidx)
##############################################	
				#P2P interact   
##############################################
dir="/shares/compbio/Group-Yang/y.wu/interaction/P2E/res/";setwd(dir)
dat=read.table("Blood_m2msmr_pass_heidi",head=T,sep="\t")
probeP=read.table("/shares/compbio/Group-Yang/y.wu/interaction/data/promoter_mprobes.list",head=F)
indx=match(dat$Outco_ID,probeP$V1,nomatch=0)
datp2p=dat[which(indx!=0),]
geneset=cbind(as.character(datp2p$Expo_Gene),as.character(datp2p$Outco_Gene))
geneset=unique(geneset)
geneset=geneset[which(geneset[,1]!=geneset[,2]),]
###############################################
     #Overlap with DNAm-gene associaitons 
###############################################
m2g=read.table("/ibscratch/users/uqywu16/omics/m2gsmr/combined_m2gsmr_pass_relaxed_heidi.smr",head=T)
m2m=read.table("/ibscratch/users/uqywu16/cSMR/P2E/combined_m2msmr_pass_heidi.smr",head=T,sep="\t")
pairsmg=paste(m2g$Expo_ID,m2g$Outco_Gene,sep="-")
pairsmm=paste(m2m$Outco_ID,m2m$Expo_Gene,sep="-")
idx=match(toupper(pairsmm),toupper(pairsmg),nomatch=0)
###############################################
     		#regional plot# 
###############################################
source("/shares/compbio/Group-Yang/y.wu/omics/regional_plot/DNAm_interact.r")
pdf("/shares/compbio/Group-Yang/y.wu/interaction/plot/new/SORT1_region.pdf",height=6,width=10)
DNAm_inter(chr=19,12875000,12915000,epi_plot=T)
DNAm_inter(chr=7,1830000,2150000,epi_plot=T)
DNAm_inter(chr=7,1830000,2280000,epi_plot=T)
DNAm_inter(chr=1,109740000,110090000,epi_plot=T)
DNAm_inter(chr=1,58606018,58199001,epi_plot=T)
DNAm_inter(chr=8,11315000,11420000,epi_plot=T)
##############################################	
			#use old resutls   
##############################################
library(data.table)
datall=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobs/ruleout/combined_m2msmr_all.smr",data.table=F)
idx = which(datall$p_SMR<(0.05/dim(datall)[1]))
smrdat = datall[idx,]
hidx = which(smrdat$p_HET>0.01)
write.table(datall[idx,],"/shares/compbio/Group-Yang/y.wu/interaction/res/blood_m2msmr_pass_smr.smr",row=F,sep="\t",qu=F)
write.table(smrdat[hidx,],"/shares/compbio/Group-Yang/y.wu/interaction/res/blood_m2msmr_pass_heidi.smr",row=F,sep="\t",qu=F)

