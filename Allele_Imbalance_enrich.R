############################################################
#				# Allele imbalance data #				   #
############################################################
# ftp://ftp.genboree.org/allelic-epigenome/json/AllelicEpigenome-sigOnly-AllDocs.json.gz
# Download the allele imbalance data from Onuchic et al. 2018
library("rjson")
filename = "/shares/compbio/Group-Yang/y.wu/interaction/AI/AllelicEpigenome-sigOnly-AllDocs.json"
result <- fromJSON(file = filename)
json_data_frame <- as.data.frame(result)
######			read json files
library(jsonlite)
result = jsonlite::fromJSON(filename)
############################################################
#	 # Extract top mQTLs (including non-significant)#	   #
############################################################
smr="/shares/compbio/Group-Yang/uqfzhan7/bin/smr_linux"
dir="/shares/compbio/Group-Yang/y.wu/interaction/AI/topmQTL/"
methy="/shares/compbio/PCTG/meQTL/LBC_BSGS/meta_meqtl_chr"
for((j=1;j<=22;j++))
do
cmd="$smr --beqtl-summary $methy${j} --peqtl-cis 0.9999 --descriptive-cis --out $dir'methy_chr'${j}'_cis_top' > $dir'methy_chr'${j}'_cis_top.log' 2>&1 "
rsub $cmd @top %50 T1 W1
done

awk 'NR==1 || FNR>1' methy_chr*_cis_top.cis.summary.txt > methy_all_cis_top.cis.summary.txt
############################################################
#			 # Allele imbalance enrich #				   #
############################################################
library(data.table)
all=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt",head=T,data.table=F)
sig=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt",head=T,sep="\t")
top=fread("/shares/compbio/Group-Yang/y.wu/interaction/AI/topmQTL/methy_all_cis_top.cis.summary.txt",head=T)
sigmprb = unique(c(as.character(sig$Expo_ID),as.character(sig$Outco_ID)))
allmprb = unique(c(as.character(all$Expo_ID),as.character(all$Outco_ID)))
sigSNPs = top$TopSNP[match(sigmprb,top$Probe,nomatch=0)]
allSNPs = top$TopSNP[match(allmprb,top$Probe,nomatch=0)]
out="/shares/compbio/Group-Yang/y.wu/interaction/AI/"
asmSNPs = as.character(read.table(paste0(out,"extract_snplist.txt"),head=F)$V3)
sigSNPs = unique(as.character(sigSNPs))
allSNPs = unique(as.character(allSNPs))

nullnum = c()
for(i in 1:1000) {
	ridx = sample(1:length(allSNPs),length(sigSNPs))
	ranSNPs = allSNPs[ridx]
	ranidx = match(ranSNPs, asmSNPs, nomatch=0)
	nullnum = c(nullnum,length(which(ranidx!=0)))
}
sigidx = match(sigSNPs, asmSNPs, nomatch=0)
obsvnum = length(which(sigidx!=0))

############################################################
#					# enrichment plot #				       #
############################################################
h=hist(nullnum,breaks=50,col=c(rep("peru",215)),xlim=c(100,250),xlab="Number of overlaps with allele-specific SNPs for DNA methylations",main="",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
abline(v=obsvnum,lty=2,lwd=2,col="red");
