
#depth: in region [0, 0.3]
#cnv: in region [0.32,0.38]
#sv: in region [0.4,0.7]


bin_depth_for_one_sample<-function(fn, startt,endd,binsize,samplename,ofn) {
  starts<-seq(from=startt,to=endd,by=binsize)
  ends<-starts+binsize
  dt<-read.table(fn,header=TRUE,sep="\t")
  ret<-rep(0,length(starts))
  for(i in 1:length(starts)) {
    rown<-which(dt[,2]>=starts[i] & dt[,2]<ends[i])
	if(length(rown>0)) {
	  ret[i]<-as.integer(median(dt[rown,7]))
	}
  }
  val<-paste(c(samplename,ret),collapse="\t")
  write.table(val,ofn,row.names=F,col.names=F,sep="\t",quote=F)  
}

prep_data_from_many_samples<-function(fns.lst,startt,endd,binsize,queue,ofnn) {
  check_jobs<-function(ooList) {
    successes<-0
    for(i in 1:length(ooList)) {
      tmp<-grep("Successfully",readLines(ooList[i]),value=TRUE)
	  if(length(tmp)==1) { successes<-successes+1 }	
    }
    return(successes)
  }

  fns<-readLines(fns.lst)
  starts<-seq(from=startt,to=endd,by=binsize)
  ends<-starts+binsize
  centers<-starts+as.integer(0.5*binsize)
  header<-paste(c("Sample",centers),collapse="\t")
  merged<-header  
  ooList<-NULL
  eoList<-NULL
  ofns<-rep("",length(fns))
  for(i in 1:length(fns)) {
    fn<-fns[i]
    samplename<-unlist(strsplit(basename(fn),split=".",fixed=TRUE))[1]
	ooList<-c(ooList,paste(c("oo.",samplename,".oo"),collapse=""))
	eoList<-c(eoList,paste(c("eo.",samplename,".eo"),collapse=""))
  }
  cmd<-paste(c("rm ", eoList, ooList),collapse=" ")
  system(cmd)
  for(i in 1:length(fns)) {
    fn<-fns[i]
	samplename<-unlist(strsplit(basename(fn),split=".",fixed=TRUE))[1]
	ofns[i]<-paste(c("bindepth_",samplename,".txt"),collapse="")	
	job<-paste(c("bsub -M8000 -q ", queue, " -eo eo.", samplename, ".eo -oo oo.",samplename,".oo Rscript plot_arch.r bindepth ", fn, " ", startt, " ", endd, " ", binsize, " ", samplename, " ", ofns[i]),collapse="")
	ooList<-c(ooList,paste(c("oo.",samplename,".oo"),collapse=""))
	eoList<-c(eoList,paste(c("eo.",samplename,".eo"),collapse=""))
	system(job)	
  }
  Sys.sleep(30)
  successes <- 0
  N<-0
  while(successes<length(ooList)) {
    print(paste(c("checking job ", N),collapse=" "))
	N<-N+1
    successes<-check_jobs(ooList)
	Sys.sleep(15)
  }
  for(i in 1:length(ofns)) {
    tmp<-readLines(ofns[i])
	merged<-c(merged,tmp)
	cmd<-paste(c("rm ", ofns[i], eoList[i], ooList[i]),collapse=" ")
	system(cmd)
  }
  write.table(merged,ofnn,row.names=F,col.names=F,sep="\t",quote=F)
  merged<-read.table(ofnn,header=TRUE,sep="\t",colClass="character",row.names=1,check.names=F)
  Samples<-row.names(merged)
  Pos<-c("Sample",colnames(merged))
  merged<-t(merged)
  merged<-rbind(Sample=Samples,merged)
  merged<-cbind(Pos,merged)
  write.table(merged,ofnn,row.names=F,col.names=FALSE,sep="\t",quote=F)
}


plot_arch_function<-function(ifn,cnvsv.ifn, ofn, refgene.ifn, IsonlyPlotDefaultGene, focusGene, padding) {
  only_default<-IsonlyPlotDefaultGene
  dt<-read.table(ifn,header=TRUE,sep="\t",colClass="numeric")
  cnvsv<-NULL
  if(cnvsv.ifn!="NULL") {
    cnvsv<-read.table(cnvsv.ifn,header=F,sep="\t",colClass="character")
  }

  xMM<-max(dt[,1])
  xmm<-min(dt[,1])
  x_len<-xMM-xmm
  my_xmm<-xmm-0.1*x_len
  my_xMM<-xMM+0.01*x_len

  refgene<-read.table(refgene.ifn,header=TRUE,sep="\t")
  focus_S<-as.numeric(refgene[refgene[,"gName"]==focusGene & refgene[,"isDefault"]=="Yes","txS"])
  focus_E<-as.numeric(refgene[refgene[,"gName"]==focusGene & refgene[,"isDefault"]=="Yes","txE"])
  mychr<-unique(refgene[refgene[,"gName"]==focusGene & refgene[,"isDefault"]=="Yes","chr"])
  if(padding>0) {
    my_xmm<-focus_S-padding
    my_xMM<-focus_E+padding
  }
  x_len<-my_xMM-my_xmm
  dt<-dt[dt[,1]>=my_xmm & dt[,1]<=my_xMM,]
  pos<-as.numeric(dt[,1])
  nSmp<-ncol(dt)-1

  overlap_interval<-(as.numeric(refgene[,"txS"])<=my_xmm & as.numeric(refgene[,"txE"])>=my_xmm) | #left boundary in
                    (as.numeric(refgene[,"txS"])<=my_xMM & as.numeric(refgene[,"txE"])>=my_xMM) | #right boundary in
                    (as.numeric(refgene[,"txS"])>=my_xmm & as.numeric(refgene[,"txE"])<=my_xMM) # gene in
  samechr<-refgene[,"chr"] %in% c(mychr)
  refgene<-refgene[samechr & overlap_interval,]  
  if(only_default==TRUE) { refgene<-refgene[refgene[,"isDefault"]=="Yes",]; }
  nIso<-nrow(refgene)

  margin<-7
  pdf(ofn,width=10,height=5+0.9*nSmp,useDingbats=FALSE)  
  par(mai = c(0.1, 0.1, 0.1, 0.1))
  plot(-100,xlim=c(my_xmm-0.1*x_len,my_xMM+0.01*x_len),ylim=c(0,nSmp+4-margin-10.5),bty="n",ylab="",xlab="",cex.axis=0.5,xaxt="n",yaxt="n")
  
  lines(c(my_xmm,my_xMM),c(-0.5-margin,-0.5-margin))
  tick_mm<-as.integer(my_xmm/100000)
  tick_MM<-as.integer(my_xMM/100000)
  ticks<-seq(from=tick_mm,to=tick_MM,by=1)/10
  ppos<-ticks*1000000
  ticks<-ticks[ppos>my_xmm & ppos<my_xMM]
  ppos<-ppos[ppos>my_xmm & ppos<my_xMM]
  for(i in 1:length(ppos)) {
    lines(c(ppos[i],ppos[i]),c(-0.5-margin,-0.7-margin))
    text(ppos[i],-0.7-margin,ticks[i],cex=0.5,pos=1)
  }
  
  text(my_xmm+0.45*x_len,nSmp-margin+2.2,mychr,cex=1,pos=1)
  lines(c(focus_S,focus_S),c(0,nSmp-margin+2),lty=2)
  lines(c(focus_E,focus_E),c(0,nSmp-margin+2),lty=2)
  YYY<-nSmp-margin+2
  yyy<-nSmp-margin+0.5
  stepp<-1.4/nIso
  yposs<-seq(from=yyy,to=YYY,by=stepp)
  for(j in 1:nrow(refgene)) {
    gName<-refgene[j,"gName"]
    isoform<-refgene[j,"nmID"]
    ort<-refgene[j,"strand"]
    txS<-as.numeric(refgene[j,"txS"])
    txE<-as.numeric(refgene[j,"txE"])
    cdsS<-as.numeric(refgene[j,"cdsS"])
    cdsE<-as.numeric(refgene[j,"cdsE"])
    exS<-as.numeric(unlist(strsplit(refgene[j,"exonS"],split=",")))
    exE<-as.numeric(unlist(strsplit(refgene[j,"exonE"],split=",")))
    lines(c(txS,txE),c(yposs[j],yposs[j]))
    for(k in 1:length(exS)) {
      myS<-exS[k]
      myE<-exE[k]
      rect(xleft=myS,xright=myE,ytop=yposs[j]+0.05*stepp,ybottom=yposs[j]-0.05*stepp)
      if(cdsS>myS & cdsS<myE) { myS<-cdsS; }
      if(cdsE>myS & cdsE<myE) { myE<-cdsE; }
      if(myS>=cdsS & myE<=cdsE) {
          rect(xleft=myS,xright=myE,ytop=yposs[j]+0.1*stepp,ybottom=yposs[j]-0.1*stepp)
      }
    }
    if(txS-my_xmm>0.2*x_len) {
      text(txS,yposs[j],paste(c(gName,isoform, ort),collapse=" | "),cex=0.5,pos=2)
    } else if(my_xMM-txE>0.2*x_len) {
      text(txE,yposs[j],paste(c(gName,isoform, ort),collapse=" | "),cex=0.5,pos=4)
    } else {
      text(my_xmm,yposs[j],paste(c(gName,isoform, ort),collapse=" | "),cex=0.5,pos=2)
    }
  } 
  
  for(i in 2:ncol(dt)) {  
    sjid<-colnames(dt)[i]
    truSJID<-unlist(strsplit(sjid,split="_"))[1]
    val_max<-max(dt[,i])  
    my_val<-(i-2-margin)+0.3*dt[,i]/val_max  # region [0, 0.3]
    lines(c(my_xmm,my_xmm),c(i-2-margin,i-2-margin+0.8))
    lines(c(my_xmm,my_xMM),c(i-2-margin,i-2-margin))
    #lines(pos,my_val)
    n_pos<-c(pos[1],pos,pos[length(pos)])
    range(n_pos)
    n_val<-c(i-2-margin,my_val,i-2-margin)
    polygon(n_pos,n_val,col="gray",lwd=0.1,border=NA)
    lines(c(my_xmm,my_xmm-0.007*x_len),c(i-2-margin+0.3,i-2-margin+0.3))
    text(my_xmm-0.002*x_len,i-2-margin+0.3,sprintf("%.1f",val_max),pos=2,cex=0.5)
    text(my_xmm,i-2-margin+0.15,sjid,cex=0.5,pos=2)
      
    if(!is.null(cnvsv)) {
      if(sjid %in% cnvsv[,1] | truSJID %in% cnvsv[,1]) {
        rown<-which(cnvsv[,1] %in% c(sjid,truSJID))

        for(j in 1:length(rown)) {
          type<-cnvsv[rown[j],2]
          marker<-cnvsv[rown[j],3]
          colorr<-cnvsv[rown[j],4]
          if(type=="CNV") {
            chr<-unlist(strsplit(marker,split=".",fixed=TRUE))[1]
            startt<-as.numeric(unlist(strsplit(marker,split=".",fixed=TRUE))[2])
            endd<-as.numeric(unlist(strsplit(marker,split=".",fixed=TRUE))[3])
            if(chr==mychr) {
              yy<-i-2-margin+0.32+runif(1,min=0,max=0.08)
              lines(c(startt,endd),c(yy,yy),col=colorr,lwd=2)
            }
          }
          if(type=="SV") {
            chrA<-unlist(strsplit(marker,split=".",fixed=TRUE))[1]
            posA<-as.numeric(unlist(strsplit(marker,split=".",fixed=TRUE))[2])
            ortA<-unlist(strsplit(marker,split=".",fixed=TRUE))[3]
            chrB<-unlist(strsplit(marker,split=".",fixed=TRUE))[4]
            posB<-as.numeric(unlist(strsplit(marker,split=".",fixed=TRUE))[5])
            ortB<-unlist(strsplit(marker,split=".",fixed=TRUE))[6]
            if(chrA==mychr && chrB==mychr) {
              yy<-i-2-margin+0.42
              yM<-yy+0.36
              middle<-0.5*posA+0.5*posB
              lines(c(posA,middle),c(yy,yM),col=colorr,lwd=2)
              lines(c(middle,posB),c(yM,yy),col=colorr,lwd=2)
              points(c(posA,posB),c(yy,yy),pch=16,cex=1,col=colorr)
            }
            if(chrA==mychr && chrB!=mychr) {
              yy<-i-2-margin+0.42
              yM<-yy+0.36            
              lines(c(posA,posA),c(yy,yM),col=colorr,lwd=2)            
              points(posA,yy,pch=16,cex=1,col=colorr)
            }
            if(chrA!=mychr && chrB==mychr) {
              yy<-i-2-margin+0.42
              yM<-yy+0.36            
              lines(c(posB,posB),c(yy,yM),col=colorr,lwd=2)            
              points(posB,yy,pch=16,cex=1,col=colorr)            
            }        
          }
        }
      }
    } 
  }
  dev.off()    
}

plot_arch_omin<-function(fnslst, cnvsv.ifn, startt, endd, binsize, refgene.ifn,focusGene, pdfofn, padding, IsonlyPlotDefaultGene,queue) {
  ofnn<-paste(c("merged_bindepth_",focusGene,".txt"),collapse="")
  prep_data_from_many_samples(fnslst, startt,endd, binsize,queue,ofnn)
  plot_arch_function(ofnn,cnvsv.ifn, pdfofn, refgene.ifn, IsonlyPlotDefaultGene, focusGene, padding)
}

#ifn<-"input.txt"
#cnvsv.ifn<-"CNVSV.txt"
#ofn<-"test.pdf"
#refgene.ifn<-"ncbiRefSeqCurated_hg19_default.txt"
#IsonlyPlotDefaultGene<-TRUE
#focusGene<-"MYCN"
#padding<-500000

need_help<-FALSE
args<-commandArgs(trailingOnly=TRUE)

if(length(args)==12) {
  if(args[1]=="omni") {
    fnslst<-args[2]
	cnvsv.ifn<-args[3]
	startt<-as.integer(args[4])
	endd<-as.integer(args[5])
	binsize<-as.integer(args[6])
	refgene.ifn<-args[7]
	focusGene<-args[8]
	pdfofn<-args[9]
	padding<-as.integer(args[10])
	IsonlyPlotDefaultGene<-args[11]
	queue<-args[12]
    plot_arch_omin(fnslst, cnvsv.ifn, startt, endd, binsize, refgene.ifn, focusGene, pdfofn, padding, IsonlyPlotDefaultGene,queue)
  } else {
    print("Usage: Rscript ~ omni <fns.lst> <cnvsv.ifn|NULL> <startt> <endd> <binsize> <refgene.ifn> <focusGene|MYCN> <pdf.ofn> <padding|-1|5000000> <IsonlyPlotDefaultGene|Y|N> <queue|compbio>")
  }
} else if(length(args)==8) {
  if(args[1]=="plotarch") {
    ifn<-args[2]
    cnvsv.ifn<-args[3]
    ofn<-args[4]
    refgene.ifn<-args[5]
    IsonlyPlotDefaultGene<-FALSE
    if(args[6] %in% c("Y", "y","yes","Yes")) {
      IsonlyPlotDefaultGene<-TRUE
    }
    focusGene<-args[7]
    padding<-as.integer(args[8])
    plot_arch_function(ifn,cnvsv.ifn, ofn, refgene.ifn, IsonlyPlotDefaultGene, focusGene, padding)
  } else {
    print("Usage: Rscript ~ plotarch <data.mat.ifn> <cnvsv.ifn|NULL> <pdf.ofn> <refgene.def.ifn> <onlyPlotDefaultIsoform|Y|N> <focusGene|MYCN> <padding|-1|5000000>")
    print("data.mat.ifn: a data matrix, row is chromosomal bins, col1 = chr, col2=pos, col3- are sample names")
    print("cnvsv.ifn: cnv/sv of samples, one record per row; three columns: samplename, type (CNV/SV), marker(chr1.2345.5678 or chr1.2345.+.chr2.4567.-.ACTG)")
    print("refgene.def.ifn: /research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/mRNA/ZL/reference/ncbiRefSeqCurated_hg*default.txt")
  }
} else if(length(args)==7) {
  if(args[1]=="bindepth") {
    ifn<-args[2]
	startt<-as.integer(args[3])
	endd<-as.integer(args[4])
	binsize<-as.integer(args[5])
	samplename<-args[6]
	ofn<-args[7]
    bin_depth_for_one_sample(ifn, startt,endd,binsize,samplename,ofn)	
  } else if(args[1]=="depthbinall") {
    fns.lst<-args[2]
	startt<-as.integer(args[3])
	endd<-as.integer(args[4])
	binsize<-as.integer(args[5])
	queue<-args[6]
	ofnn<-args[7]    
    prep_data_from_many_samples(fns.lst, startt,endd, binsize,queue,ofnn)
  } else {
    print("Usage: Rscript ~ depthbinall <fns.lst> <startt> <endd> <binsize> <queue> <ofn>")
    print("fns.lst: sequencerr output file list; each file has 7 columns: chr pos A C G T N")
    print("startt: start position to visualize; endd: end position; binsize: size of bin")    
  } 
} else {
  print("Usage: Rscript ~ plotarch|depthbinall\n")
  print("Usage: Rscript ~ plotarch <data.mat.ifn> <cnvsv.ifn|NULL> <pdf.ofn> <refgene.def.ifn> <onlyPlotDefaultIsoform|Y|N> <focusGene|MYCN> <padding|-1|5000000>\n")  
  print("Usage: Rscript ~ depthbinall <fns.lst> <startt> <endd> <binsize> <queue> <ofn>\n")
  print("Usage: Rscript ~ omni <fns.lst> <cnvsv.ifn|NULL> <startt> <endd> <binsize> <refgene.ifn> <focusGene|MYCN> <pdf.ofn> <padding|-1|5000000> <IsonlyPlotDefaultGene|Y|N> <queue|compbio>")
} 
