# Concordance among HiC loops from rao et al., Jin et al.and Jung et al.#
######################################################################### 
        # replicate Hi-C among three datasets (all sig. loops) #
#########################################################################
pgltools="/shares/compbio/Group-Yang/y.wu/interaction/HiC/pgltools-master/sh/pgltools"
junghic="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_sigloops"
raohic="/shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_primary"
jinhic="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC"
jungraorep="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_rao_rep"
jungjinrep="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_jin_rep"
jinraorep="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/jin_rao_rep"
#############   two pairs
$pgltools intersect -a $junghic.pgl -b $raohic.pgl -wo > $jungraorep.pgld
$pgltools intersect -a $junghic.pgl -b $jinhic.pgl -wo > $jungjinrep.pgld
$pgltools intersect -a $jinhic.pgl -b $raohic.pgl -wo > $jinraorep.pgld
#############   three pairs
jungovlapjin="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_ovlap_jin"
jinovlapjung="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/jin_ovlap_Jung"
raoovlapjin="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/rao_ovlap_jin"
jungovlapjinovlaprao="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_ovlap_jin_ovlap_rao"
jinovlapjungovlaprao="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/jin_ovlap_Jung_ovlap_rao"
raoovlapjinovlapjung="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/rao_ovlap_jin_ovlap_jung"
awk '{print $1,$2,$3,$4,$5,$6}' $jungjinrep.pgld | sort | uniq > $jungovlapjin.bedpe
awk '{print $7,$8,$9,$10,$11,$12}' $jungjinrep.pgld | sort | uniq > $jinovlapjung.bedpe
awk '{print $7,$8,$9,$10,$11,$12}' $jinraorep.pgld | sort | uniq > $raoovlapjin.bedpe
$pgltools formatbedpe $jungovlapjin.bedpe > $jungovlapjin.pgl
$pgltools formatbedpe $jinovlapjung.bedpe > $jinovlapjung.pgl
$pgltools formatbedpe $raoovlapjin.bedpe > $raoovlapjin.pgl
$pgltools intersect -a $jungovlapjin.pgl -b $raohic.pgl -wo > $jungovlapjinovlaprao.pgld
$pgltools intersect -a $jinovlapjung.pgl -b $raohic.pgl -wo > $jinovlapjungovlaprao.pgld
$pgltools intersect -a $raoovlapjin.pgl -b $junghic.pgl -wo > $raoovlapjinovlapjung.pgld
############	summarize
raoall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_primary.pgl | wc -l",intern=T))
jinall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC.pgl | wc -l",intern=T))
jungall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_sigloops.pgl | wc -l",intern=T))
raorep3 = as.numeric(system("awk '{print $1,$2,$3,$4,$5,$6}' /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/rao_ovlap_jin_ovlap_jung.pgld | sort | uniq | wc -l",intern=T))
jinrep3 = as.numeric(system("awk '{print $1,$2,$3,$4,$5,$6}' /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/jin_ovlap_Jung_ovlap_rao.pgld | sort | uniq | wc -l",intern=T))
jungrep3 = as.numeric(system("awk '{print $1,$2,$3,$4,$5,$6}' /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_ovlap_jin_ovlap_rao.pgld | sort | uniq | wc -l",intern=T))
round(raorep3/raoall,4)
round(jinrep3/jinall,4)
round(jungrep3/jungall,4)
#########################################################################


######################################################################### 
         # replicate Hi-C among three datasets (top 10K loops) #
#########################################################################
#############   top 10K
Jung = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_sigloops.pgl",data.table=F)
Jin = fread("/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC.pgl",data.table=F)
Jungidx = which(Jung$V7 >= sort(Jung$V7,decreasing=T)[9448])
Jinidx = which(Jin$V7 <= sort(Jin$V7)[9448])
write.table(Jung[Jungidx,],"/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_top94K_loops.pgl",row=F,col=F,sep="\t",qu=F)
write.table(Jin[Jinidx,],"/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC_top94K_loops.pgl",row=F,col=F,sep="\t",qu=F)
############# 	overlap
pgltools="/shares/compbio/Group-Yang/y.wu/interaction/HiC/pgltools-master/sh/pgltools"
junghic="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_top94K_loops"
raohic="/shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_10Kloops"
jinhic="/shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC_top94K_loops"
jungraorep="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_rao_rep"
jungjinrep="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_jin_rep"
jinraorep="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/jin_rao_rep"
#############   two pairs
$pgltools intersect -a $junghic.pgl -b $raohic.pgl -wo > $jungraorep.pgld
$pgltools intersect -a $junghic.pgl -b $jinhic.pgl -wo > $jungjinrep.pgld
$pgltools intersect -a $jinhic.pgl -b $raohic.pgl -wo > $jinraorep.pgld
#############   three pairs
jungovlapjin="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_ovlap_jin"
jinovlapjung="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/jin_ovlap_Jung"
raoovlapjin="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/rao_ovlap_jin"
jungovlapjinovlaprao="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_ovlap_jin_ovlap_rao"
jinovlapjungovlaprao="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/jin_ovlap_Jung_ovlap_rao"
raoovlapjinovlapjung="/shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/rao_ovlap_jin_ovlap_jung"
awk '{print $1,$2,$3,$4,$5,$6}' $jungjinrep.pgld | sort | uniq > $jungovlapjin.bedpe
awk '{print $7,$8,$9,$10,$11,$12}' $jungjinrep.pgld | sort | uniq > $jinovlapjung.bedpe
awk '{print $7,$8,$9,$10,$11,$12}' $jinraorep.pgld | sort | uniq > $raoovlapjin.bedpe
$pgltools formatbedpe $jungovlapjin.bedpe > $jungovlapjin.pgl
$pgltools formatbedpe $jinovlapjung.bedpe > $jinovlapjung.pgl
$pgltools formatbedpe $raoovlapjin.bedpe > $raoovlapjin.pgl
$pgltools intersect -a $jungovlapjin.pgl -b $raohic.pgl -wo > $jungovlapjinovlaprao.pgld
$pgltools intersect -a $jinovlapjung.pgl -b $raohic.pgl -wo > $jinovlapjungovlaprao.pgld
$pgltools intersect -a $raoovlapjin.pgl -b $junghic.pgl -wo > $raoovlapjinovlapjung.pgld
############	summarize
raoall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_10Kloops.pgl | wc -l",intern=T))
jinall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC_top94K_loops.pgl | wc -l",intern=T))
jungall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_top94K_loops.pgl | wc -l",intern=T))
raorep3 = as.numeric(system("awk '{print $1,$2,$3,$4,$5,$6}' /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/rao_ovlap_jin_ovlap_jung.pgld | sort | uniq | wc -l",intern=T))
jinrep3 = as.numeric(system("awk '{print $1,$2,$3,$4,$5,$6}' /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/jin_ovlap_Jung_ovlap_rao.pgld | sort | uniq | wc -l",intern=T))
jungrep3 = as.numeric(system("awk '{print $1,$2,$3,$4,$5,$6}' /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_ovlap_jin_ovlap_rao.pgld | sort | uniq | wc -l",intern=T))
round(raorep3/raoall,4)
round(jinrep3/jinall,4)
round(jungrep3/jungall,4)


#########################################################################
			 # Jaccard Index plot #
#########################################################################
############	all significant loops from each study
raoall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_primary.pgl | wc -l",intern=T))
jinall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC.pgl | wc -l",intern=T))
jungall = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_sigloops.pgl | wc -l",intern=T))
Jung_rao_ol1 = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_rao_rep.pgld | wc -l",intern=T))
Jung_jin_ol1 = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/Jung_jin_rep.pgld | wc -l",intern=T))
jin_rao_ol1 = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/jin_rao_rep.pgld | wc -l",intern=T))
Jung_rao_JI1 = Jung_rao_ol1/(jungall+raoall-Jung_rao_ol1)
Jung_jin_JI1 = Jung_jin_ol1/(jungall+jinall-Jung_jin_ol1)
jin_rao_JI1 = jin_rao_ol1/(jinall+raoall-jin_rao_ol1)
###########   top 10K loops from each study
raotop = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/RaoData/GM12878_10Kloops.pgl | wc -l",intern=T))
jintop = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JinData/JinHiC_top94K_loops.pgl | wc -l",intern=T))
jungtop = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/JungData/Jung_PCHiC_GM_top94K_loops.pgl | wc -l",intern=T))
Jung_rao_ol2 = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_rao_rep.pgld | wc -l",intern=T))
Jung_jin_ol2 = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/Jung_jin_rep.pgld | wc -l",intern=T))
jin_rao_ol2 = as.numeric(system("cat /shares/compbio/Group-Yang/y.wu/interaction/HiC/hicrep/hicrep1/jin_rao_rep.pgld | wc -l",intern=T))
Jung_rao_JI2 = Jung_rao_ol2/(jungtop+raotop-Jung_rao_ol2)
Jung_jin_JI2 = Jung_jin_ol2/(jungtop+jintop-Jung_jin_ol2)
jin_rao_JI2 = jin_rao_ol2/(jintop+raotop-jin_rao_ol2)
###########	barplot
pdf("/shares/compbio/Group-Yang/y.wu/interaction/HiC/HiC_replicates_datasets_JI.pdf")
JIall = c(Jung_rao_JI1,Jung_jin_JI1,jin_rao_JI1,Jung_rao_JI2,Jung_jin_JI2,jin_rao_JI2)
plot<-barplot(height=JIall,col=c("#8DD3C7","darkolivegreen4","navy","orange","purple","#377EB8"),axes=TRUE, axisnames=F,ylim=c(0,0.3),cex.lab=1.5,cex.axis=1.5,names.arg=names,cex.names=1.5,ylab="Jaccard Index",xlab="")
legend("topleft",c(expression(paste("Rao et al. vs. Jung et al. (",italic(P), " < 0.01)")),expression(paste("Jung et al. (",italic(P), " < 0.01) vs. Jin et al. (FDR < 0.1)")),"Rao et al. vs. Jin et al. (FDR < 0.1)","Rao et al. vs. Jung et al. (High-confidence)","Jin et al. (High-confidence) vs. Jung et al. (High-confidence)","Rao et al. vs. Jin et al. (High-confidence)"),pt.cex=1,cex=1,fill=c("#8DD3C7","darkolivegreen4","navy","orange","purple","#377EB8"),bty="n")
dev.off()
#########################################################################
			      # end #
#########################################################################
