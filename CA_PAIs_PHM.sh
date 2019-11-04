################################################################################## 
				# configuration #				   			     
##################################################################################
SMR="/shares/compbio/Group-Yang/uqfzhan7/bin/smr_linux"
MLTSMR="/shares/compbio/Group-Yang/y.wu/program/mtsmr/mtSMR0.10/mltsmr"
K9AC="/shares/compbio/Group-Yang/t.qi/data/ROSMAP/besd/haQTL_all"
K27AC="/shares/compbio/Group-Yang/y.wu/blueprint/besd/tcel_K27AC_log2rpm_peer_10_all_summary"
K4ME1="/shares/compbio/Group-Yang/y.wu/blueprint/besd/tcel_K4ME1_log2rpm_peer_10_all_summary"
HQTL="/shares/compbio/Group-Yang/y.wu/database/besd/hQTL_merged"
CAQTL="/shares/compbio/Group-Yang/y.wu/database/besd/caQTL"
METHY="/shares/compbio/PCTG/meQTL/LBC_BSGS/meta_meqtl_chr"
CAGE="/shares/compbio/Group-Yang/eQTL/CAGE/CAGE.cis.chr"
eQTLGen="/shares/compbio/Group-Yang/y.wu/SMRdata/eQTLGen/cis-eQTLsFDR0.05-ProbeLevel.txt_besd"
PQTL="/shares/compbio/Group-Yang/y.wu/pQTL/besdall/pQTL_3KProt_sparse"
GWASDIRT='/shares/compbio/Group-Yang/y.wu/SMRdata/public_meta'
PHENO=("" "height" "bmi" "whradjbmi" "hdl" "ldl" "tg" "ea" "ra" "scz" "cad" "t2d" "cd" "uc" "ad")
GWASDATA=(""  "Yengo_height.raw" "Yengo_BMI.raw" "whradjbmi_giant_2015.raw" "hdl_glgc_2013.raw " "ldl_glgc_2013.raw" "tg_glgc_2013.raw" "GWAS_EA3_excl23andMe.raw" "ra_okada_2014.raw" "scz_pgc_2014.raw" "cad_cardiogram_2015.raw" "meta_DGU_whole_sum_for_SMR.txt" "cd_iibdgc_2015.raw" "uc_iibdgc_2015.raw" "alzheimers_igap_2013.raw")
NPHENO=$(($(echo ${#PHENO[@]}) - 1 ))
REFERENCE="/shares/compbio/Group-Yang/reference_data/hrs/hrs_1kg_hwe1e-6"
INDI="/shares/compbio/Group-Yang/uqzzhu1/dataset/SepIDs/hrs_hm3_m01_hwe_1e-6_u05_geno.grm.id"
HRSSNPLIST="/shares/compbio/Group-Yang/t.qi/data/SepIDs/hrs_1kg_impRsq_0.3.snplist"
ICLPROBES="/shares/compbio/Group-Yang/y.wu/heidi/heidi2/result/LBS_DNAm_mqtl.list"
ICLOUTCO="/shares/compbio/Group-Yang/y.wu/heidi/heidi2/result/molecular_outcomes_with_sigQtl.list"
EXPROBES="/shares/compbio/Group-Yang/y.wu/SMRdata/MHC_probes.list"
GENELIST="/shares/compbio/Group-Yang/y.wu/database/refseq_isoform_SMR.txt"
DIR="/shares/compbio/Group-Yang/y.wu/allSMR/"
######################################################## 
			 # C2C #						   
########################################################
OUTPUT=$DIR"c2csmr/"
for((i=1;i<=22;i++))
do
rsub "${SMR} --beqtl-summary ${CAQTL} --beqtl-summary ${CAQTL} \
 --bfile ${REFERENCE}_chr${i} --extract-snp ${HRSSNPLIST} --keep ${INDI}\
 --maf 0.01 --diff-freq-prop 0.5 --cis-wind 2000 --out ${OUTPUT}c2csmr_chr${i} --thread-num 3 \
 > ${OUTPUT}c2csmr_chr${i}.log 2>&1" @c2csmr %8 T3 W168
done
############ lower the threshold 1e-3
OUTPUT=$DIR"c2csmr/c2ctest/"
for((i=1;i<=22;i++))
do
rsub "${SMR} --beqtl-summary ${CAQTL} --beqtl-summary ${CAQTL} --peqtl-smr 1e-3 \
 --bfile ${REFERENCE}_chr${i} --extract-snp ${HRSSNPLIST} --keep ${INDI}\
 --maf 0.01 --diff-freq-prop 0.5 --cis-wind 500 --out ${OUTPUT}c2csmr_chr${i} --thread-num 3 \
 > ${OUTPUT}c2csmr_chr${i}.log 2>&1" @c2csmr %8 T3 W168
done
########################################################
	    	  # PAIs replicated in PHM #
########################################################
library(data.table)
c2c = read.table("/shares/compbio/Group-Yang/y.wu/allSMR/c2csmr/c2csmr_pass_heidi.txt",head=T,sep="\t")
PP = fread("gunzip -c /shares/compbio/Group-Yang/y.wu/database/caQTL/posterior_prob.tsv.gz",head=T)
sigidx = which(PP$Causality_j_2_k>0.5 | PP$Causality_k_2_j>0.5 | PP$Pleiotropy>0.5)
PPsig = PP[sigidx,]
pair1 = paste(c2c$Expo_ID,c2c$Outco_ID,sep="-")
pair2 = paste(c2c$Outco_ID,c2c$Expo_ID,sep="-")
pairall = paste(PP$Peak_j,PP$Peak_k,sep="-")
pairsig = paste(PPsig$Peak_j,PPsig$Peak_k,sep="-")
aidx1 = match(pair1,pairall,nomatch=0)
aidx2 = match(pair2,pairall,nomatch=0)
midx1 = match(pair1,pairsig,nomatch=0)
midx2 = match(pair2,pairsig,nomatch=0)
aidx = unique(c(aidx1, aidx2))
midx = unique(c(midx1, midx2))
length(which(midx!=0))/length(which(aidx!=0))
########################################################
	    	 # PHM replicated in SMR&HEIDI #
########################################################
library(data.table)
c2c = read.table("/shares/compbio/Group-Yang/y.wu/allSMR/c2csmr/c2csmr_all.txt",head=T,sep="\t")
PP = fread("gunzip -c /shares/compbio/Group-Yang/y.wu/database/caQTL/posterior_prob.tsv.gz",head=T)
sigidx = which(PP$Causality_j_2_k>0.5 | PP$Causality_k_2_j>0.5)
PPsig = PP[sigidx,]
pairc2c = paste(c2c$Expo_ID,c2c$Outco_ID,sep="-")
pair1 = paste(PPsig$Peak_k,PPsig$Peak_j,sep="-")
pair2 = paste(PPsig$Peak_j,PPsig$Peak_k,sep="-")
midx1 = match(pair1,pairc2c,nomatch=0)
midx2 = match(pair2,pairc2c,nomatch=0)
midx = which(midx1==0 & midx2!=0)
midx1[midx] = midx2[midx]
c2crepall = c2c[midx1,]
sigidx = which(c2crepall$p_SMR <= 0.05/nrow(c2crepall) & c2crepall$p_HEIDI >= 0.01)
length(sigidx)/nrow(c2crepall)
########################################################
