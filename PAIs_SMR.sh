################################################################# 
			# DNAm interactions
#################################################################
COHORT="Blood"
SMR="/shares/compbio/Group-Yang/uqfzhan7/bin/smr_linux"
METHY="/shares/compbio/PCTG/meQTL/LBC_BSGS/meta_meqtl_chr"
ProProbes="/shares/compbio/Group-Yang/y.wu/allSMR/m2msmr/P2E/promoter_mprobes.list"
REFERENCE="/shares/compbio/Group-Yang/reference_data/hrs/hrs_1kg_hwe1e-6"
INDI="/shares/compbio/Group-Yang/uqzzhu1/dataset/SepIDs/hrs_hm3_m01_hwe_1e-6_u05_geno.grm.id"
HRSSNPLIST="/shares/compbio/Group-Yang/t.qi/data/SepIDs/hrs_1kg_impRsq_0.3.snplist"
EXPROBES="/shares/compbio/Group-Yang/y.wu/omic2/QC/excl_mhc_tqtl_mprbs.list"
EXSNPS="/shares/compbio/Group-Yang/y.wu/SMRdata/MHC_snps.list"
GENELIST="/shares/compbio/Group-Yang/y.wu/database/refseq_isoform_SMR.txt"
#################################################################  
			 # P2E (including)  
#################################################################
OUTPUT="/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/"
for((i=1;i<=22;i++))
do
rsub "${SMR} --beqtl-summary ${METHY}${i} --extract-exposure-probe ${ProProbes} \
 --beqtl-summary ${METHY}${i} \
 --bfile ${REFERENCE}_chr${i} --exclude-snp ${EXSNPS} --extract-snp ${HRSSNPLIST} --keep ${INDI} \
 --maf 0.01 --cis-wind 2000 --out ${OUTPUT}${COHORT}_m2msmr_chr${i} --thread-num 20 \
 > ${OUTPUT}${COHORT}_m2msmr_chr${i}.log 2>&1" @interaction %4 T7 W150
done
################################################################# 
     	      # Exclude DNAm pairs in same promoter 
#################################################################
promp=read.table("/shares/compbio/Group-Yang/y.wu/interaction/data/methylation_promoters.txt",head=F)
megeo=read.table("/shares/compbio/Group-Yang/y.wu/database/GPL16304-47833.txt",skip=22,head=T,sep="\t")
indir="/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/"
outdir="/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/"
for (i in 1:22) {
	print(i)
	dat=read.delim(paste(indir,"Blood_m2msmr_chr",i,".smr",sep=""),head=T,sep="\t")
	indx=match(dat$Expo_ID,promp$V2,nomatch=0)
	tmp=which(dat$Outco_bp>promp$V5[indx] & dat$Outco_bp<=promp$V6[indx])
	data=dat[-tmp,]
	expoindx=match(data$Expo_ID,megeo$ID,nomatch=0)
	outcindx=match(data$Outco_ID,megeo$ID,nomatch=0)
	data$Expo_Gene=megeo$Closest_TSS_gene_name[expoindx]
	data$Outco_Gene=megeo$Closest_TSS_gene_name[outcindx]
	write.table(data,paste(outdir,"Blood_m2msmr_chr",i,".smr",sep=""),row=F,qu=F,sep="\t")
}
#################################################################
			   # summarize  
#################################################################
setwd("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/")
datall=c()
for(i in 1:22) {
	print(i)
	file=paste0("Blood_m2msmr_chr",i,".smr")
	if(file.exists(file)==1) {
	data=read.table(file,head=T,sep="\t",colClass=c(rep(c("character","numeric"),5),"numeric",rep("character",2),rep("numeric",12)))
	datall=rbind(datall,data)
	}
}
write.table(datall,"/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all.txt",row=F,sep="\t",qu=F)
#################################################################
		          # remove MHC  
#################################################################
library(data.table)
all = fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all.txt",head=T,data.table=F,sep="\t",stringsAsFactors=F)
#####			MHC region
mhcStart = 26056121; mhcEnd = 33377699
allmhc = which( (all$Expo_Chr==6 & all$Expo_bp>=mhcStart & all$Expo_bp<=mhcEnd) | (all$Outco_Chr==6 & all$Outco_bp>=mhcStart & all$Outco_bp<=mhcEnd) )
all = all[-allmhc,]
outfile2 = "/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc.txt"
write.table(all,outfile2,row=F,sep="\t",qu=F)
################################################################ 
   		    # remove duplicated pairs  
################################################################
library(data.table)
outdir = "/shares/compbio/Group-Yang/y.wu/interaction/final/"
all = fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc.txt",head=T,data.table=F)
allpair1 = paste(all$Expo_ID, all$Outco_ID, sep="-")
allpair2 = paste(all$Outco_ID, all$Expo_ID, sep="-")
alldupidx = match(allpair2,allpair1,nomatch=0); 

allrmidx = c(); index1=alldupidx[which(alldupidx!=0)]; index2=which(alldupidx!=0)
for(i in 1:length(index1)) {
	idx1 = index1[i]; 	idx2 = index2[i]
	if((all[idx1,]$b_Expo/all[idx1,]$se_Expo)^2 >= (all[idx2,]$b_Expo/all[idx2,]$se_Expo)^2) 
	{
		allrmidx = c(allrmidx,idx2)
	}
}
all = all[-allrmidx,]
sigidx = which(all$p_SMR <=0.05/nrow(all) & all$p_HEIDI >=0.01)
sig = all[sigidx,]
outfile1 = paste0(outdir,"Blood_m2msmr_all_exclmhc_rmdup.txt")
outfile2 = paste0(outdir,"Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt")
write.table(all,outfile1,row=F,sep="\t",qu=F)
write.table(sig,outfile2,row=F,sep="\t",qu=F)
################################################################# 
			    # end #  
#################################################################
