####################################################################
   # APA enrichment for predicted loops from different methods #
####################################################################
#################################		
   		# Data format #
#################################		
library(data.table); options(scipen=999)
setwd("/shares/compbio/Group-Yang/y.wu/interaction/APA/tools/")
loops = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_10Kloops.pgl",,head=F,data.table=F)
sigDNAm = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/sigDNAm_001.pgl",head=F,data.table=F)
corrDNAm = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/sample/sigCorr.pgl",head=F,data.table=F)
CAcorr = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/sample/sigCACorr.pgl",head=F,data.table=F)
sigCASMR = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/sample/sigCASMR.pgl",head=F,data.table=F)
sigCAPHM = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/sample/sigCAPHM.pgl",head=F,data.table=F)
ranDNAm = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/Partition/ransigDNAm.pgl",head=F,data.table=F)
loops = loops[,c(1:6)]; loops$color="255,0,0"; loops$comment = "null"
sigDNAm = sigDNAm[,c(1:6)]; sigDNAm$color="255,0,0"; sigDNAm$comment = "null"
corrDNAm = corrDNAm[,c(1:6)]; corrDNAm$color="255,0,0"; corrDNAm$comment = "null"
CAcorr = CAcorr[,c(1:6)]; CAcorr$color="255,0,0"; CAcorr$comment = "null"
sigCASMR = sigCASMR[,c(1:6)]; sigCASMR$color="255,0,0"; sigCASMR$comment = "null"
sigCAPHM = sigCAPHM[,c(1:6)]; sigCAPHM$color="255,0,0"; sigCAPHM$comment = "null"
ranDNAm = ranDNAm[,c(1:6)]; ranDNAm$color="255,0,0"; ranDNAm$comment = "null"
colnames(loops)[1:6] = colnames(sigDNAm)[1:6] = colnames(corrDNAm)[1:6] = colnames(CAcorr)[1:6] = 
colnames(sigCASMR)[1:6] = colnames(sigCAPHM)[1:6] = colnames(ranDNAm)[1:6] = c("chr1","x1","x2","chr2","y1","y2")
write.table(unique(loops),"Rao_loops.txt",row=F,sep="\t",qu=F)
write.table(unique(sigDNAm),"sigDNAm.txt",row=F,sep="\t",qu=F)
write.table(unique(corrDNAm),"corrDNAm.txt",row=F,sep="\t",qu=F)
write.table(unique(CAcorr),"CAcorr.txt",row=F,sep="\t",qu=F)
write.table(unique(sigCASMR),"sigCASMR.txt",row=F,sep="\t",qu=F)
write.table(unique(sigCAPHM),"sigCAPHM.txt",row=F,sep="\t",qu=F)
write.table(unique(ranDNAm),"ranDNAm.txt",row=F,sep="\t",qu=F)
################################		
   		# APA enrich #
################################		
######	 input
DIR="/shares/compbio/Group-Yang/y.wu/interaction/APA/tools/"
tools=$DIR"juicer_tools_1.14.08.jar"
sigDNAm=$DIR"sigDNAm.txt"
corrDNAm=$DIR"corrDNAm.txt"
CAcorr=$DIR"CAcorr.txt"
sigCASMR=$DIR"sigCASMR.txt"
sigCAPHM=$DIR"sigCAPHM.txt"
ranDNAm=$DIR"ranDNAm.txt"
######	 output
ressigDNAm=$DIR"res_sigDNAm"
rescorrDNAm=$DIR"res_corrDNAm"
resCAcorr=$DIR"res_CAcorr"
ressigCASMR=$DIR"res_sigCASMR"
ressigCAPHM=$DIR"res_sigCAPHM"
resranDNAm=$DIR"res_ranDNAm"
######	 analysis
Rscript ~/bin/rsub.R java -jar $tools apa -u -r 5000 https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic $sigDNAm $ressigDNAm @sigDNAm %60 T1 W48
Rscript ~/bin/rsub.R java -jar $tools apa -u -r 5000 https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic $corrDNAm $rescorrDNAm @corrDNAm %60 T1 W48
Rscript ~/bin/rsub.R java -jar $tools apa -u -r 5000 https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic $CAcorr $resCAcorr @CAcorr %60 T1 W48
Rscript ~/bin/rsub.R java -jar $tools apa -u -r 5000 https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic $sigCASMR $ressigCASMR @sigCASMR %60 T1 W48
Rscript ~/bin/rsub.R java -jar $tools apa -u -r 5000 https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic $sigCAPHM $ressigCAPHM @sigCAPHM %60 T1 W48
Rscript ~/bin/rsub.R java -jar $tools apa -u -r 5000 https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic $ranDNAm $resranDNAm @ranDNAm %60 T1 W48
###############################		
   	# Summarize measures#
###############################
DIR="/shares/compbio/Group-Yang/y.wu/interaction/APA/tools/"
sigDNAm_APA=read.table(paste0(DIR,"res_sigDNAm","/5000/gw/measures.txt"),head=F)$V2[4]
corrDNAm_APA=read.table(paste0(DIR,"res_corrDNAm","/5000/gw/measures.txt"),head=F)$V2[4]
CAcorr_APA=read.table(paste0(DIR,"res_CAcorr","/5000/gw/measures.txt"),head=F)$V2[4]
sigCASMR_APA=read.table(paste0(DIR,"res_sigCASMR","/5000/gw/measures.txt"),head=F)$V2[4]
sigCAPHM_APA=read.table(paste0(DIR,"res_sigCAPHM","/5000/gw/measures.txt"),head=F)$V2[4]
sigCAPHMC_APA=read.table(paste0(DIR,"res_sigCAPHMC","/5000/gw/measures.txt"),head=F)$V2[4]
ranDNAm_APA=read.table(paste0(DIR,"res_ranDNAm","/5000/gw/measures.txt"),head=F)$V2[4]

sigDNAm_Zscore=read.table(paste0(DIR,"res_sigDNAm","/5000/gw/measures.txt"),head=F)$V2[6]
corrDNAm_Zscore=read.table(paste0(DIR,"res_corrDNAm","/5000/gw/measures.txt"),head=F)$V2[6]
CAcorr_Zscore=read.table(paste0(DIR,"res_CAcorr","/5000/gw/measures.txt"),head=F)$V2[6]
sigCASMR_Zscore=read.table(paste0(DIR,"res_sigCASMR","/5000/gw/measures.txt"),head=F)$V2[6]
sigCAPHM_Zscore=read.table(paste0(DIR,"res_sigCAPHM","/5000/gw/measures.txt"),head=F)$V2[6]
ranDNAm_Zscore=read.table(paste0(DIR,"res_ranDNAm","/5000/gw/measures.txt"),head=F)$V2[6]
####################################################################
   							# End #
####################################################################