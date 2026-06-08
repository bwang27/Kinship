#03/19/2025 Xiaotu Ma
trimed.median<-function(xx) {
  cutt<-0
  dnn<-quantile(xx,0.01)
  upp<-quantile(xx,0.99)
  if(sum(xx>dnn & xx<upp)>2) {
    use<-(xx>dnn & xx<upp)
    cutt<-quantile(xx[use],0.95)
  }
  return(as.numeric(cutt))
}
beutify.it.0<-function(xx,offsett=1e6) {
  zz<-ceiling((10*offsett+xx)/10^floor(log10(xx)))
  return(zz*(10^floor(log10(xx))))
}
plot.one.panel<-function(poss,MAF,CVG,medCVG,cnvsv,xoffset,chrlen, nroww, upoffset.low,upoffset.high,nPrev,spacer,shrinkage=1,seq.type="WGS",chr,figmod) {
  intervals<-floor(chrlen/1e7)
  if(!is.null(CVG)) {
    CVG[CVG>medCVG]<-medCVG
    if(length(poss)>0) { 
      for(j in 1:length(poss)) {
        lines(c(xoffset+poss[j],xoffset+poss[j]),c(upoffset.high,upoffset.high+shrinkage*CVG[j]/medCVG),col="orange")
      }
    }
  }
  cex.points<-0.4
  if(seq.type=="WES") cex.points<-1
  if(!is.null(MAF)) {
    MAF[MAF>=1]<-1
    MAF[MAF<=0]<-0
    points(xoffset+poss,upoffset.high+shrinkage*MAF,pch=16,cex=cex.points,col="#555555")
  }
  the.col<-"skyblue"
  if(nPrev%%2==1) the.col<-"springgreen"  
  rect(xleft=xoffset,xright=xoffset+chrlen,ytop=upoffset.high,ybottom=upoffset.low,col=the.col,border=NA)
  if(nPrev>0) {
    rect(xleft=xoffset-spacer*0.7,xright=xoffset-spacer*0.3,ytop=upoffset.high+shrinkage,ybottom=upoffset.low,col="red",border=NA)
  }
  for(j in 0:intervals) {
    lines(c(xoffset+j*1e7,xoffset+j*1e7),c(upoffset.low,upoffset.high),lwd=4)
    tick.len<-upoffset.high-upoffset.low
    for(k in 1:9) {
      if((j+k/10)*1e7<chrlen) {
        tick.height<-0.8
        if(k%%2==0) {
          tick.height<-0.5
        }
        lines(c(xoffset+(j+k/10)*1e7,xoffset+(j+k/10)*1e7),c(upoffset.low+tick.len*tick.height,upoffset.high),lwd=4)
      }
    }
  }
  if(!is.null(cnvsv)) {
    if(nrow(cnvsv)>=1) {
      for(j in 1:nrow(cnvsv)) {
        theCol<-"red"
        if("Color" %in% colnames(cnvsv)) {
          theCol<-as.character(cnvsv[j,"Color"])
        }
        if(cnvsv[j,"Type"] %in% c("CNV","CNVM","LOH")) {
          theChr<-unlist(strsplit(cnvsv[j,"Marker"],split=".",fixed=TRUE))[1]
          if(theChr==chr) {
            theStart<-as.integer(unlist(strsplit(cnvsv[j,"Marker"],split=".",fixed=TRUE))[2])
            theEnd<-as.integer(unlist(strsplit(cnvsv[j,"Marker"],split=".",fixed=TRUE))[3])
            if(figmod %in% c("D","N")) {
              shifts<-runif(n=1,max=0.45,min=-0.45)
              lines(c(xoffset+theStart,xoffset+theEnd),c(upoffset.high+shrinkage*(0.5+shifts),upoffset.high+shrinkage*(0.5+shifts)),lwd=8,col=theCol)
            }
            if(figmod %in% c("DN")) {
               shifts<-runif(n=1,max=0.15,min=-0.15)
               lines(c(xoffset+theStart,xoffset+theEnd),c(upoffset.high+shrinkage*(0.8+shifts),upoffset.high+shrinkage*(0.8+shifts)),lwd=8,col=theCol)
            }
          }
        }
        if(cnvsv[j,"Type"] %in% c("SV")) {      
          vals<-unlist(strsplit(cnvsv[j,"Marker"],split=".",fixed=TRUE))
          if(vals[1]==chr) {
            thePos<-as.integer(vals[2])
            if(figmod %in% c("D","N")) {
              shifts<-runif(n=1,max=0.45,min=-0.45)        
              points(xoffset+thePos,upoffset.high+shrinkage*(0.5+shifts),pch=1,col=theCol,cex=4,lwd=4)
            }
            if(figmod %in% c("DN")) {
              shifts<-runif(n=1,max=0.15,min=-0.5)
              points(xoffset+thePos,upoffset.high+shrinkage*(0.8+shifts),pch=1,col=theCol,cex=4,lwd=4)
            }
          }
          if(vals[4]==chr) {
            thePos<-as.integer(vals[5])
            if(figmod %in% c("D","N")) {
              shifts<-runif(n=1,max=0.4,min=-0.4)        
              points(xoffset+thePos,upoffset.high+shrinkage*(0.5+shifts),pch=1,col=theCol,cex=4,lwd=4)
            }
            if(figmod %in% c("DN")) {
              shifts<-runif(n=1,max=0.15,min=-0.5)
              points(xoffset+thePos,upoffset.high+shrinkage*(0.8+shifts),pch=1,col=theCol,cex=4,lwd=4)
            }
          }
        }      
      }
    }
  }
  xpos<-xoffset+(1:intervals)*1e7
  return(xpos)
}

arrange.panels<-function(chr.size.ifn,dt0,figmod,cnvsv,maxcvg) {
  chr.sz<-read.table(chr.size.ifn,header=T,sep="\t",colClass=c("character",rep("numeric",3)))
  nSNP<-nrow(dt0)
  seq.type<-"WGS"
  if(nSNP<6e5) {
    seq.type<-"WES"
  }
  chrs<-apply(cbind("chr",c(1:22,"X","Y")),1,paste,collapse="")
  rows<-sort(unique(chr.sz[,"Row"]))
  max.len<-0
  panel.height<-1
  panel.dist<-0.1
  letter.track.height<-0.8
  track.height<-3*(panel.height+panel.dist)+letter.track.height
  spacer<-3e6
  shrinkage<-0.9
  for(i in 1:length(rows)) {
    len<-sum(chr.sz[chr.sz[,"Row"]==rows[i],"Length"])
    if(len>max.len) max.len<-len
  }  
  nroww<-length(rows)*track.height
  maxY<-nroww
  maxX<-max.len*1.05
  max_med <- 0
  if(figmod %in% c("D","D1","DN","D1N1")) {
    max_med<-trimed.median(dt0[,"TinD"])*1.5
  }
  if(figmod %in% c("N","N1","DN","D1N1")) {
    medd.G<-trimed.median(dt0[,"TinN"])*1.5
    if(medd.G>max_med) max_med<-medd.G
  }
  if(maxcvg!=-1) {
    max_med<-maxcvg
  }
  if(max_med<5) max_med<-5
  #medd<-beutify.it.0(max_med,offsett=0)
  medd<-ceiling(max_med)
  plot(-5,xlim=c(0,maxX),ylim=c(0,maxY),xlab="",ylab="",main="",xaxt="n",yaxt="n")
  for(i in 1:length(chrs)) {
    the.row<-chr.sz[chr.sz[,1]==chrs[i],"Row"]
    the.Y<-(length(rows)-the.row)*track.height
    starts<-0
    nPrev<-0
    nFut<-0
    if(i>1) {
      prev_chrs<-chrs[1:(i-1)]
      prev_chrs<-intersect(prev_chrs, chr.sz[chr.sz[,"Row"]==chr.sz[i,"Row"],"Chr"])
      nPrev<-length(prev_chrs)
      starts<-sum(chr.sz[as.character(chr.sz[,1]) %in% prev_chrs,"Length"])
      if(length(prev_chrs)>0)   starts<-starts+length(prev_chrs)*spacer
      fut_chrs<-setdiff(chrs[i:length(chrs)],chrs[i])
      fut_chrs<-intersect(fut_chrs, chr.sz[chr.sz[,"Row"]==chr.sz[i,"Row"],"Chr"])
      nFut<-length(fut_chrs)
      if(nFut==0) {
        themy.x<-starts+chr.sz[chr.sz[,"Chr"]==chrs[i],"Length"] + spacer
        panel.lbl<-c("D-G","G","D")
        yy<-(length(rows)-the.row)*track.height+letter.track.height+(0:2)*(panel.height+panel.dist)+panel.dist
        for(kk in 1:length(yy)) {
          lines(c(themy.x, themy.x), c(yy[kk],yy[kk]+shrinkage),lwd=3)
          for(ll in 0:5) {
            lines(c(themy.x,themy.x+maxX*0.002),c(yy[kk]+shrinkage*ll/5,yy[kk]+shrinkage*ll/5),lwd=3)
            text(themy.x+maxX*0.002,yy[kk]+shrinkage*ll/5,medd*ll/5,cex=3,pos=4)
          }
          text(themy.x+maxX*0.01,yy[kk]+shrinkage*2.5/5,panel.lbl[kk],cex=4,pos=4)
        }
      }
    }
    poss0<-dt0[as.character(dt0[,"Chr"])==chrs[i],"Pos"]
    ord<-order(as.numeric(poss0))
    poss<-poss0[ord]
    if(figmod %in% c("DN","D1N1")) {
      dCVG<-dt0[as.character(dt0[,"Chr"])==chrs[i],"TinD"][ord]
      dMAF<-dt0[as.character(dt0[,"Chr"])==chrs[i],"dMAF"][ord]
      nMAF<-dt0[as.character(dt0[,"Chr"])==chrs[i],"nMAF"][ord]
      nCVG<-dt0[as.character(dt0[,"Chr"])==chrs[i],"TinN"][ord]
      diffMAF<-abs(dMAF-nMAF)  
      plot.one.panel(poss=poss, MAF=dMAF,CVG=dCVG,medCVG=medd,cnvsv=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+2*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+2*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
      plot.one.panel(poss=poss, MAF=nMAF,CVG=nCVG,medCVG=medd,cnvsv=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+1*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+1*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
      xpos<-plot.one.panel(poss=poss, MAF=diffMAF,CVG=NULL,medCVG=medd,cnvsv=cnvsv,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+0*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+0*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
    }
  if(figmod %in% c("D","D1")) {
      dCVG<-dt0[as.character(dt0[,"Chr"])==chrs[i],"TinD"][ord]
      dMAF<-dt0[as.character(dt0[,"Chr"])==chrs[i],"dMAF"][ord]
      plot.one.panel(poss=poss, MAF=NULL,CVG=dCVG,medCVG=medd,cnvsv=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+2*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+2*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
      plot.one.panel(poss=poss, MAF=dMAF,CVG=NULL,medCVG=medd,cnvsv=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+1*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+1*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
      xpos<-plot.one.panel(poss=poss, MAF=NULL,CVG=NULL,medCVG=medd,cnvsv=cnvsv,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+0*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+0*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
    }
  if(figmod %in% c("N","N1")) {
      nMAF<-dt0[as.character(dt0[,"Chr"])==chrs[i],"nMAF"][ord]
      nCVG<-dt0[as.character(dt0[,"Chr"])==chrs[i],"TinN"][ord]  
      plot.one.panel(poss=poss, MAF=NULL,CVG=nCVG,medCVG=medd,cnvsv=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+2*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+2*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
      plot.one.panel(poss=poss, MAF=nMAF,CVG=NULL,medCVG=medd,cnvsv=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+1*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+1*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
      xpos<-plot.one.panel(poss=poss, MAF=NULL,CVG=NULL,medCVG=medd,cnvsv=cnvsv,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+0*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+0*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type,chr=chrs[i],figmod)
    }
    for(j in 1:length(xpos)) {
      text(xpos[j], the.Y+letter.track.height-0.1, labels=j*10,cex=3)
    }
    mid<-median(xpos)
    text(mid, the.Y+letter.track.height*0.5,chrs[i],cex=7)
  }
  for(j in 1:length(rows)) {
    the.x<-0-5e6
    yy<-(j-1)*track.height+letter.track.height+(0:2)*(panel.height+panel.dist)+panel.dist
    for(k in 1:length(yy)) {
      lines(c(the.x, the.x), c(yy[k],yy[k]+shrinkage),lwd=3)
      for(ll in 0:5) {
        lines(c(the.x-maxX*0.002,the.x),c(yy[k]+shrinkage*ll/5,yy[k]+shrinkage*ll/5),lwd=3)
        text(the.x-maxX*0.002,yy[k]+shrinkage*ll/5,ll/5,cex=3,pos=2)
      }
    }
  }
}

read.cnvsv<-function(cnvsv.ifn) {
  cnvsv<-NULL
  if(!is.null(cnvsv.ifn) & cnvsv.ifn!="NULL") {
    cnvsv<-read.table(cnvsv.ifn,header=T,colClass="character",sep="\t")
    colnames(cnvsv)[1:3]<-c("Sample","Type","Marker")
    if(ncol(cnvsv)>=4) colnames(cnvsv)[4]<-"Color"  
    for(i in 1:nrow(cnvsv)) {
      cnvsv[i,"Marker"]<-gsub("chr23","chrX",cnvsv[i,"Marker"])
      cnvsv[i,"Marker"]<-gsub("chr24","chrY",cnvsv[i,"Marker"])
    }
  }  
  return(cnvsv)
}

do.one.file<-function(chr.inf.ifn, dt.ifn, figmod,cnvsv.ifn,jpg.ofn,sample.name,maxcvg) {
  minNdepth<-15
  cnvsv<-read.cnvsv(cnvsv.ifn)
  if(!is.null(cnvsv)) {
    if(sample.name %in% cnvsv[,"Sample"]) {
      cnvsv<-cnvsv[cnvsv[,"Sample"]==sample.name,]
    } else {
      cnvsv<-NULL
    }
  }
  if(figmod %in% c("D","DN","D1","N","N1")) {
    dt0<-read.table(dt.ifn,header=T,nrow=0-1)
    nms<-colnames(dt0)
    if(figmod=="D1") {
      nms[nms=="nM"]<-"MinD"
      nms[nms=="nT"]<-"TinD"
    }
    if(figmod=="N1") {
      nms[nms=="nM"]<-"MinN"
      nms[nms=="nT"]<-"TinN"
    }
    colnames(dt0)<-nms
  }
  if(figmod=="D1N1") {
    fns<-unlist(strsplit(dt.ifn,split="=",fixed=TRUE))
    if(length(fns)==2) {
      dt.ifn_D<-fns[1]
      dt.ifn_N<-fns[2]
      dt_D<-read.table(dt.ifn_D,header=TRUE,sep="\t")
      dt_N<-read.table(dt.ifn_N,header=TRUE,sep="\t")
      if(length(intersect(c("Chr","Pos","nM","nT"),colnames(dt_D)))<4 | length(intersect(c("Chr","Pos","nM","nT"),colnames(dt_N)))<4) {
        print("the D1 file and N1 file need to have column names Chr, Pos, nM, and nT")
        return;
      }
      row.names(dt_D)<-apply(cbind(dt_D[,"Chr"],as.character(dt_D[,"Pos"])),1,paste,collapse=".")
      row.names(dt_N)<-apply(cbind(dt_N[,"Chr"],as.character(dt_N[,"Pos"])),1,paste,collapse=".")
      common<-intersect(row.names(dt_D),row.names(dt_N))
      dt_D<-dt_D[common,]
      dt_N<-dt_N[common,]
      nmsD<-colnames(dt_D)
      nmsD[nmsD=="nM"]<-"MinD"
      nmsD[nmsD=="nT"]<-"TinD"
      colnames(dt_D)<-nmsD
      nmsN<-colnames(dt_N)
      nmsN[nmsN=="nM"]<-"MinN"
      nmsN[nmsN=="nT"]<-"TinN"
      colnames(dt_N)<-nmsN
      dt0<-cbind(dt_D,dt_N[,c("MinN","TinN")])
    } else {
      print("please input two data.ifns in the format of ifnD=ifnN")
      return;
    }
  }
  if(length(intersect(c("MinD","TinD","dMAF"),colnames(dt0)))<2 & figmod %in% c("D","DN","D1","D1N1")) {
    print("Input file needs columns TinD and MinD|dMAF");
    return;  
  }
  if(length(intersect(c("MinN","TinN","nMAF"),colnames(dt0)))<2 & figmod %in% c("N","DN","N1","D1N1")) {
    print("Input file needs columns TinN and MinN|nMAF");
    return;  
  }
  if(figmod %in% c("D","DN","D1N1","D1")) {
    if(!("dMAF" %in% colnames(dt0))) {
      dMAF<-dt0[,"MinD"]/(0.001+dt0[,"TinD"])
      dt0<-cbind(dt0,dMAF)
    }
  }
  if(figmod %in% c("N","DN","D1N1","N1")) {
    if(!("nMAF" %in% colnames(dt0))) {
      nMAF<-dt0[,"MinN"]/(0.001+dt0[,"TinN"])
      dt0<-cbind(dt0,nMAF)
    }
  }
  if("TinN" %in% colnames(dt0)) {
    xx<-dt0$TinN
    if(maxcvg==-1) {
      if(sum(xx>minNdepth)>100) {
        upp<-as.numeric(quantile(xx[xx>10],0.98))
        cutofff<-mean(xx[xx<upp & xx>10]) + 4*sd(xx[xx<upp & xx>10])  
        dt0<-dt0[dt0$TinN>minNdepth & dt0$TinN<cutofff,]
      }
      if(sum(xx>1)<5) {
        dt0<-dt0[dt0$TinD>40,]
      }
  }
    if(maxcvg>1) {
      dt0<-dt0[dt0$TinN>10,]
    }
  }
  jpeg(jpg.ofn,width=8000,height=4000)
  arrange.panels(chr.size.ifn=chr.inf.ifn,dt0,figmod,cnvsv,maxcvg)
  dev.off()
}

args<-commandArgs(TRUE)
#print(args)
if(length(args)==6) {
  do.one.file(chr.inf.ifn=args[1], dt.ifn=args[2], figmod=args[3], cnvsv.ifn=args[4], jpg.ofn=args[5],sample.name=args[6],maxcvg=-1)
}

if(length(args)==7) {
  do.one.file(chr.inf.ifn=args[1], dt.ifn=args[2], figmod=args[3], cnvsv.ifn=args[4], jpg.ofn=args[5],sample.name=args[6],maxcvg=as.integer(args[7]))
}


if(length(args)!=6 & length(args)!=7) {
  print("Usage: Rscript ~ <chr.size.ifn> <data.ifn|=N1> <figmod|D|N|DN|D1|N1|D1N1> <cnvsv.ifn|NULL> <jpg.ofn> <sample.name> | <max_cvg>")
  print("cnvsv file: <Sample> <Type(CNV|SV)> <Marker(chr.start.end or chrA.pA.+.chrB.pB.+.)> <Color>")
  print("data.ifn: can be DG.ifn or D1.ifn=G1.ifn when figmog=D1N1")
  print("data.ifn: must have Chr Pos | MinD TinD | MinN | TinN | nM nT columns")
}

