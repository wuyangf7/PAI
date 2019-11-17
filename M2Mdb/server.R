library(shiny)
#library(OmicCircos)
me2tr=readRDS("metro.rds");
m2m=readRDS("m2m.rds");
#glist=readRDS("db/glist_hg19.rds");
#meIndx=readRDS("db/meIndx.rds");
#eIndx=readRDS("db/eIndx.rds");
#meSMR=readRDS("db/me2eSMR.rds")
plotWind=1000000;
# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x);

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x);

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x);

shinyServer(function(input, output) {
  
  output$mymapping = DT::renderDataTable(DT::datatable({ 
    #data <- me2tr[,c(1,2,26,4:6,27,8,23:24)]
    data <- m2m[,c(1,2,3,4,5,6,7,8,23:24)]
    if(input$Chr != "All") {
      data <- data[which(data$Expo_Chr == input$Chr),]
    }
    if(input$meprobeid != ""  ) {
      data <- data[which(trim(input$meprobeid) == data$Expo_ID),]
    }
    if(input$eprobeid != "" ) {
      data <- data[which(data$Outco_ID == trim(input$eprobeid)),]
    }
    if(input$geneid != "" ) {
      data <- data[which(data$Expo_Gene == trim(input$geneid)),]
    }
    if(input$psmr != "" ) {
      data <- data[which(data$p_SMR <= as.numeric(trim(input$psmr))),]
    }
    if(input$phet != "" ) {
      data <- data[which(data$p_HET >= as.numeric(trim(input$phet))),]
    }
    names(data)=c("Exposure DNAm","Exposure DNAm chr","Exposure DNAm gene","Exposure DNAm BP","Outcome DNAm","Outcome DNAm chr","Outcome DNAm gene","Outcome DNAm BP","pSMR","pHEIDI")
    curChr=unique(data[,1])
    data$pHEIDI=format(data$pHEIDI,  digits = 3,  scientific=TRUE)
    data$pSMR=format(data$pSMR,  digits = 3,  scientific=TRUE) 
    data
    
  })) 
  output$text1 = renderUI({
    if(trim(input$meprobeid) == "" | trim(input$eprobeid) == "")
      {
        paste("Please specify one mProbe and one eProbe to plot!");
      } else {
        str1=" Regional circle plot:";
        str2=paste(" the outter circle represents meQTL of probe ",input$meprobeid );
        str3=paste(" the middle circle represents eQTL of probe ",input$eprobeid,sep="");
        str4=paste(" the inner circle represents SMR of probe ",input$eprobeid," to methylation probe ",input$meprobeid," and nearby methylation probes.",sep="");
        HTML(paste(str1, str2, str3, str4,sep = '<br/>'))
      }
    
  })
  output$about = renderUI({
      str1="<b> Integrative analysis of omics data identifies putative functional genes and DNA elements for complex traits </b>";
      str2=" Yang Wu, Jian Zeng, Futao Zhang, Zhihong Zhu, Ting Qi, Zhili Zheng, Luke Lloyd-Jones, Naomi Wray, Peter M. Visscher, Allan F. McRae, Jian Yang* . ";
      str3=" <br/>If you would like to use any of the results displayed or downloaded from this page in a publication, please cite this work.";
      str4=" <br/>For further questions, please contact the corresponding authors at jian.yang@imb.uq.edu.au ";
      str5=" <br/>This app was developed by Futao Zhang.";
      HTML(paste(str1, str2, str3, str4, str5, sep = '<br/>'))
  
  })
  output$myplot = renderPlot({
    if(trim(input$meprobeid) != "" & trim(input$eprobeid) != "") {
      idx=which(eIndx[,2]==trim(input$eprobeid));
      if(length(idx)==1) {
        eprobechr=as.numeric(eIndx[idx,1]);
        eprobebp=as.numeric(eIndx[idx,3]);
        #get eQTL of this region
        ename=paste("db/eQTL_chr",eprobechr,".rds",sep="");
        eQTLdata=readRDS(ename);
        idx4=which(eIndx[,1]==eprobechr);
        subeIndx=eIndx[idx4,];
        idx3=which(subeIndx[,2]==trim(input$eprobeid));
        eqtl=eQTLdata[[idx3]];
        
        idx1=which(glist[,1]==eprobechr & glist[,2]<=eprobebp+plotWind & glist[,3]>= eprobebp-plotWind);
        #adjust gene list
        if(max(glist[idx1,3])<(eprobebp+plotWind) & glist[idx1[length(idx1)]+1,1]==eprobechr) idx1=c(idx1,idx1[length(idx1)]+1);
        if(min(glist[idx1,2])>(eprobebp-plotWind) & glist[idx1[1]-1,1]==eprobechr) idx1=c(idx1[1]-1,idx1);
        seg.f=glist[idx1,];
        seg.f=seg.f[order(seg.f[,2]),];
        # get eQTL of this region
        idx7=which(eqtl[,1]<=seg.f[dim(seg.f)[1],3] & eqtl[,1]>= seg.f[1,2])
        eqtl=eqtl[idx7,]
        pvalue_eqtl=cbind(eprobechr,eqtl[,c(1,1,2)]);
        pvalue_gene=c()
        for(i in 1:dim(pvalue_eqtl)[1])
        {
          idx2=which(pvalue_eqtl[i,2]>=seg.f[,2]) 
          pvalue_gene=c(pvalue_gene,as.character(seg.f[idx2[length(idx2)],4]))
        }
        pvalue_eqtl[,1]=pvalue_gene;
        colnames(pvalue_eqtl)=c("egene","po","SNP","pvalue");
        
        #get meQTL of this region
        idx4=which(meIndx[,1]==eprobechr);
        submeIndx=meIndx[idx4,];
        idx3=which(submeIndx[,2]==trim(input$meprobeid));
        mename=paste("db/meQTL_chr",submeIndx[idx3,1],".rds",sep="");
        meQTLdata=readRDS(mename);
        meqtl=meQTLdata[[idx3]];
        idx6=which(meqtl[,1]<=seg.f[dim(seg.f)[1],3] & meqtl[,1]>= seg.f[1,2])
        meqtl=meqtl[idx6,]
        pvalue_meqtl=cbind(eprobechr,meqtl[,c(1,1,2)]);
        mgene=c()
        for(i in 1:dim(pvalue_meqtl)[1])
        {
          idx5=which(pvalue_meqtl[i,2]>=seg.f[,2]) 
          mgene=c(mgene,as.character(seg.f[idx5[length(idx5)],4]))
        }
        pvalue_meqtl[,1]=mgene;
        colnames(pvalue_meqtl)=c("mgene","po","SNP","pvalue");

        #get me2eSMR
        idx8=which(meSMR[,3]==trim(input$eprobeid));
        msmr=meSMR[idx8,];
        #idx9=which(msmr[,5]<=seg.f[dim(seg.f)[1],3] & msmr[,5]>= seg.f[1,2]) # using top SNP bp
        idx9=which(msmr[,2]<=seg.f[dim(seg.f)[1],3] & msmr[,2]>= seg.f[1,2]) # using methyl probe bp
        msmr=msmr[idx9,]
        #msmr=msmr[order(msmr[,5]),]
        msmr=msmr[order(msmr[,2]),] # using methyl probe bp
        #pvalue_me2e=cbind(eprobechr, msmr[,c(5,1,6:7)]);
        pvalue_me2e=cbind(eprobechr, msmr[,c(2,1,6:7)]); # using methyl probe bp
        geneid=c()
        for(i in 1:dim(pvalue_me2e)[1])
        {
          idxX=which(seg.f[,2]<=pvalue_me2e[i,2])
          geneid=c(geneid,as.character(seg.f[idxX[length(idxX)],4]))
        }
        pvalue_me2e[,1]=geneid;
        colnames(pvalue_me2e)=c("gene","po","probe","pvalue","pheidi")
        
        seg.f[,1]=seg.f[,4]
        seg.name=unique(as.character(seg.f[,1]));  
        seg.num=length(seg.name);
        seg.f=cbind(seg.f,1);
        db <- segAnglePo(seg.f,seg=seg.name);
        colors <- rainbow(seg.num,alpha=0.5);
        colors2 <- rainbow(10,alpha=0.8);
        par(mar=c(2,2,2,2));
        plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="");
        circos(R=400,type="chr",cir=db,print.chr.lab=TRUE, W=4,scale=TRUE,col=colors );
        circos(R=330,cir=db,W=70, mapping=pvalue_meqtl,col.v=4,type="s",B=TRUE,lwd=1,col="red"); # meQTL
        circos(R=260,cir=db,W=70, mapping=pvalue_eqtl,col.v=4,type="s",B=FALSE,lwd=1,col="blue"); #eQTL
        circos(R=190,cir=db,W=70, mapping=pvalue_me2e,col.v=4,type="s",B=FALSE,lwd=1,col="purple"); #me2eSMR
        #lableinfo=pvalue_me2e[order(pvalue_me2e[,4],decreasing=TRUE),]
        #lableinfo=lableinfo[which(lableinfo[,5]>=0.05),]
        #lableinfo=lableinfo[which(lableinfo[,3]!=trim(input$meprobeid)),]
        #lables=lableinfo[1:3,]
        #circos(R=190,cir=db , W=50, mapping=lables ,type="label",side="in" , col=c( " black " ) , cex=1 ) ;
        
        lableinfo=pvalue_me2e[which(pvalue_me2e[,3]==trim(input$meprobeid)),]
        lables=lableinfo[,1:3]
        circos(R=190,cir=db , W=50, mapping=lables ,type="label",side="in" , col=c( " blue "  ) , cex=1 ) ;
        
        
      } else print("Wrong probe name!");
   
    }
    
  })
})
