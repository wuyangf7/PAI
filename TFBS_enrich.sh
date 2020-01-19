############################################################
#					 #TFBS formating					   #
############################################################
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/
# Download the TFBS from link above and the lists in Table S3 from Rao et al. 2014, cell. 
# Format the PAIs as a 10Kb region centered by the DNAm sites
############################################################
#					 #PAI formating						   #
############################################################
DIR="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/"; setwd(DIR)
sig=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt",head=T,sep="\t")
options(scipen=999)
sig$Expo_Chr = paste0("chr",sig$Expo_Chr); sig$Outco_Chr = paste0("chr",sig$Outco_Chr)
sigexpo = data.frame(sig$Expo_Chr, abs(sig$Expo_bp - 5000), sig$Expo_bp + 5000)
sigoutco = data.frame(sig$Outco_Chr, abs(sig$Outco_bp - 5000), sig$Outco_bp + 5000)
colnames(sigexpo) = colnames(sigoutco) = paste0("V",1:3)
sigfor = unique(rbind(sigexpo,sigoutco))
write.table(sigfor,"sigDNAm.bed",row=F,,col=F,sep="\t",qu=F)
############################################################
#					 #Random sampling					   #
############################################################
library(data.table)
out="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/RsamplesV2/"
dir="/shares/compbio/Group-Yang/reference_data/1000GP_Phase3/plinkbed/"
all=fread("/shares/compbio/Group-Yang/y.wu/interaction/P2E/jobsraw/ruleout/Blood_m2msmr_all_exclmhc_rmdup.txt",head=T,data.table=F)
sig=read.table("/shares/compbio/Group-Yang/y.wu/interaction/res/Blood_m2msmr_pass_heidi_exclmhc_rmdup.txt",head=T,sep="\t")
######			all DNAm format
all$Expo_Chr=paste0("chr",all$Expo_Chr); all$Outco_Chr=paste0("chr",all$Outco_Chr)
allexpo = data.frame(all$Expo_Chr, abs(all$Expo_bp - 5000), all$Expo_bp + 5000)
alloutco = data.frame(all$Outco_Chr, abs(all$Outco_bp - 5000), all$Outco_bp + 5000)
colnames(allexpo) = colnames(alloutco) = paste0("V",1:3)
allfor = unique(rbind(allexpo,alloutco))
######		  all genomic positions
SNPs = fread(paste0(dir,"1000G_phase3_20130502_all.bim"),head=F)
SNPs$V1 = paste0("chr",SNPs$V1);
SNPs = data.frame(SNPs$V1, abs(SNPs$V4 - 5000), SNPs$V4 + 5000)
######		  random sampling
for(i in 1:1000) {
	print(i)
	######		random DNAm
	ranidx = sample(1:nrow(allfor),21787)
	datfor = allfor[ranidx,]
	write.table(datfor,paste0(out,"nullDNAm_",i,".bed"),row=F,,col=F,sep="\t",qu=F)
	######		random GP
	idx = sample(1:nrow(SNPs),21787)
	output = SNPs[idx,]
	write.table(output,paste0(out,"nullSNPs_",i,".bed"),row=F,col=F,sep="\t",qu=F)
}
############################################################
#					 # bedtools match #					   #
############################################################
bedtools="/shares/compbio/Group-Yang/t.qi/bin/bedtools"
dir="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/"
sigDNAmbed=$dir"sigDNAm.bed"
CTCFbed=$dir"CTCF_TFBS.bed"
RAD21bed=$dir"RAD21_TFBS.bed"
ZNF143bed=$dir"ZNF143_TFBS.bed"
YY1bed=$dir"YY1_TFBS.bed"
${bedtools} intersect -a ${sigDNAmbed} -b ${CTCFbed} -wa > ${dir}sigDNAm_in_CTCF.txt
${bedtools} intersect -a ${sigDNAmbed} -b ${RAD21bed} -wa > ${dir}sigDNAm_in_RAD21.txt
${bedtools} intersect -a ${sigDNAmbed} -b ${ZNF143bed} -wa > ${dir}sigDNAm_in_ZNF143.txt
${bedtools} intersect -a ${sigDNAmbed} -b ${YY1bed} -wa > ${dir}sigDNAm_in_YY1.txt
##########################  
###random sample match ###
##########################
for((i=1;i<=1000;i++))
do
echo $i
nullDNAm="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/RsamplesV2/nullDNAm_"$i
nullSNPs="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/RsamplesV2/nullSNPs_"$i
$bedtools intersect -a $nullDNAm.bed -b ${CTCFbed} -wa > $nullDNAm"_in_CTCF.txt"
$bedtools intersect -a $nullSNPs.bed -b ${CTCFbed} -wa > $nullSNPs"_in_CTCF.txt"
$bedtools intersect -a $nullDNAm.bed -b ${RAD21bed} -wa > $nullDNAm"_in_RAD21.txt"
$bedtools intersect -a $nullSNPs.bed -b ${RAD21bed} -wa > $nullSNPs"_in_RAD21.txt"
$bedtools intersect -a $nullDNAm.bed -b ${ZNF143bed} -wa > $nullDNAm"_in_ZNF143.txt"
$bedtools intersect -a $nullSNPs.bed -b ${ZNF143bed} -wa > $nullSNPs"_in_ZNF143.txt"
$bedtools intersect -a $nullDNAm.bed -b ${YY1bed} -wa > $nullDNAm"_in_YY1.txt"
$bedtools intersect -a $nullSNPs.bed -b ${YY1bed} -wa > $nullSNPs"_in_YY1.txt"
done
#############################################################
#					#summarize enrichment					#
#############################################################
DNAmCTCF="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_CTCF.txt"
DNAmRAD21="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_RAD21.txt"
DNAmZNF143="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_ZNF143.txt"
DNAmYY1="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_YY1.txt"
SNPsCTCF="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_CTCF.txt"
SNPsRAD21="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_RAD21.txt"
SNPsZNF143="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_ZNF143.txt"
SNPsYY1="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_YY1.txt"
for((i=1; i<=1000; i++))
do
nullDNAm="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/RsamplesV2/nullDNAm_"$i
nullSNPs="/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/RsamplesV2/nullSNPs_"$i
cat $nullDNAm"_in_CTCF.txt" | uniq | wc -l >> $DNAmCTCF
cat $nullDNAm"_in_RAD21.txt" | uniq | wc -l >> $DNAmRAD21
cat $nullDNAm"_in_ZNF143.txt" | uniq | wc -l >> $DNAmZNF143
cat $nullDNAm"_in_YY1.txt" | uniq | wc -l >> $DNAmYY1
cat $nullSNPs"_in_CTCF.txt" | uniq | wc -l >> $SNPsCTCF
cat $nullSNPs"_in_RAD21.txt" | uniq | wc -l >> $SNPsRAD21
cat $nullSNPs"_in_ZNF143.txt" | uniq | wc -l >> $SNPsZNF143
cat $nullSNPs"_in_YY1.txt" | uniq | wc -l >> $SNPsYY1
done
############################################################
#						# enrichment #					   #
############################################################
sigDNAmCTCF=9454; sigDNAmRAD21=7588; 
sigDNAmZNF143=6854; sigDNAmYY1=9477;
DNAmCTCF=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_CTCF.txt",head=F)$V1
DNAmRAD21=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_RAD21.txt",head=F)$V1
DNAmZNF143=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_ZNF143.txt",head=F)$V1
DNAmYY1=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullDNAm_in_YY1.txt",head=F)$V1
SNPsCTCF=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_CTCF.txt",head=F)$V1
SNPsRAD21=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_RAD21.txt",head=F)$V1
SNPsZNF143=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_ZNF143.txt",head=F)$V1
SNPsYY1=read.table("/shares/compbio/Group-Yang/y.wu/interaction/ChIPseq/resultV2/nullSNPs_in_YY1.txt",head=F)$V1
######			  quantify
CTCF_DNAm = sigDNAmCTCF/mean(DNAmCTCF); CTCF_SNPs = sigDNAmCTCF/mean(SNPsCTCF) 
RAD21_DNAm = sigDNAmRAD21/mean(DNAmRAD21); RAD21_SNPs = sigDNAmRAD21/mean(SNPsRAD21) 
ZNF143_DNAm = sigDNAmZNF143/mean(DNAmZNF143); ZNF143_SNPs = sigDNAmZNF143/mean(SNPsZNF143) 
YY1_DNAm = sigDNAmYY1/mean(DNAmYY1); YY1_SNPs = sigDNAmYY1/mean(SNPsYY1)
############################################################
#					# combined plot #				       #
############################################################
######			TFBS enrichment
data = cbind(c(CTCF_DNAm,RAD21_DNAm,ZNF143_DNAm,YY1_DNAm), c(CTCF_SNPs,RAD21_SNPs,ZNF143_SNPs,YY1_SNPs))
datase = cbind(c(sd(DNAmCTCF)/sigDNAmCTCF,sd(DNAmRAD21)/sigDNAmRAD21,sd(DNAmZNF143)/sigDNAmZNF143,sd(DNAmYY1)/sigDNAmYY1), 
	c(sd(SNPsCTCF)/sigDNAmCTCF,sd(SNPsRAD21)/sigDNAmRAD21,sd(SNPsZNF143)/sigDNAmZNF143,sd(SNPsYY1)/sigDNAmYY1))
plot<-barplot(t(data),col=c("orange","navy"),axes=TRUE, axisnames=F,ylim=c(0,6),cex.lab=1.5,beside=T,cex.names=1.5,ylab="Fold enrichment",xlab="DNA-binding proteins",cex.axis=1.5,main="",cex.main=1.5)
legend("topleft",c("Random controls from the tested DNAm sites","Random controls from the whole genome sites"),cex=1.5,fill=c("orange","navy"),bty="n")
names = c("CTCF","RAD21","ZNF143","YY1")
axis(1, at=c(2,5,8,11), labels=names, srt=45,las=1, cex.axis=1.5)
segments(plot,t(data+datase),plot,t(data-datase), lwd=1.5)
arrows(plot,t(data+datase),plot,t(data-datase), angle=90, code=3, length=0.03)
abline(h=1,lty=2,col="red")
mtext("a",side=3,line=0.5,adj=-0.1,cex=2.5,col="black",font=1)
############################################################
#						# End #							   #
############################################################