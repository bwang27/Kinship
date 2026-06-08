
add_shared<-function(mind,tind,minr,tinr,maxCvg,maxLen,anno,typee) {
  vafD<-as.numeric(mind)/(as.numeric(tind)+0.1)
  vafR<-as.numeric(minr)/(as.numeric(tinr)+0.1)
  lenD<-maxLen*tind/maxCvg
  lenR<-maxLen*tinr/maxCvg
  if(tind>maxCvg) {    lenD<-maxLen  }
  if(tinr>maxCvg) {    lenR<-maxLen  }
  if(as.character(anno)!="") {
    text(vafD,vafR,anno,pos=4,cex=0.5,col="red")
    lines(c(vafD-lenD*0.5,vafD+lenD*0.5),c(vafR,vafR),col="red")
    lines(c(vafD,vafD),c(vafR-lenR*0.5,vafR+lenR*0.5),col="red")
  
  } else {
    lines(c(vafD-lenD*0.5,vafD+lenD*0.5),c(vafR,vafR))
    lines(c(vafD,vafD),c(vafR-lenR*0.5,vafR+lenR*0.5))
  }
  if(typee=="indel") {     points(vafD,vafR,pch=1,cex=1)  }
}


add_Ronly<-function(mind, tind,minr,tinr,maxCvg,maxLen,leftX,rightX,anno,typee) {
  vafR<-as.numeric(minr)/(as.numeric(tinr)+0.1)
  theX<-tinr/maxCvg
  tinr[which(tinr >= 490 )] = 490
  tind[which(tind >= 490 )] = 490
  if(tinr>maxCvg) theX<-1
  theX<-theX*(rightX-leftX)+leftX  
  lenD<-maxLen*tind/maxCvg
  if(tind>maxCvg) {   lenD<-maxLen  }
  if(as.character(anno)!="") {
    text(theX,vafR,anno,pos=4,cex=0.5,col="red")
    lines(c(theX-0.5*lenD,theX+0.5*lenD),c(vafR,vafR),col="red")
    lines(c(theX,theX),c(vafR+0.005,vafR-0.005),col="red")
  } else {
    lines(c(theX-0.5*lenD,theX+0.5*lenD),c(vafR,vafR))
    lines(c(theX,theX),c(vafR+0.005,vafR-0.005))
  }
  if(typee=="indel") {     points(theX,vafR,pch=1,cex=1)  }
}


add_Donly<-function(mind, tind,minr,tinr,maxCvg,maxLen,bottomY,topY,anno,typee) {
  vafD<-as.numeric(mind)/(as.numeric(tind)+0.1)
  tinr[which(tinr >= 490 )] = 490
  tind[which(tind >= 490 )] = 490
  theY<-tind/maxCvg
  if(tinr>maxCvg) theY<-1
  theY<-theY*(topY-bottomY)+bottomY
  lenR<-maxLen*tinr/maxCvg
  if(tind>maxCvg) {   lenR<-maxLen  }
  if(as.character(anno)!="") {
    text(vafD,theY,anno,pos=4,cex=0.5,col="red")
    lines(c(vafD,vafD),c(theY-0.5*lenR,theY+0.5*lenR),col="red")
    lines(c(vafD+0.005,vafD-0.005),c(theY,theY),col="red")
  } else {
    lines(c(vafD,vafD),c(theY-0.5*lenR,theY+0.5*lenR))
    lines(c(vafD+0.005,vafD-0.005),c(theY,theY))
  }
  if(typee=="indel") {     points(vafD,theY,pch=1,cex=1)  }
}

maf2d_one<-function(onesample,sjid,cloudfn,info=sjid) {
  maxCvg<-500
  maxLen<-0.03
  mini_offset<-0.02
  xyoffset<-0-0.3
  plot(-100,xlim=c(xyoffset,1.05),ylim=c(xyoffset,1.05),xlab="",ylab="",main="" ,bty="n",xaxt="n",yaxt="n",las=1)
  
  if(cloudfn!="") {
    ifn<-cloudfn
	if(file.exists(ifn)) {
	  kk<-readLines(ifn)
	  nr<-length(grep("Single",kk,value=TRUE))
	  dronly<-read.table(ifn,header=F,nrow=nr,sep="\t",row.names=1)
	  shared<-read.table(ifn,header=F,skip=nr,sep="\t",row.names=1)
	  colnames(shared)<-row.names(shared)
	  colfunc <- colorRampPalette(c("#FFFFFF", "#AAAAAA"))
	  nCol<-100
	  cols<-colfunc(nCol)
	  grids<-row.names(shared)
	  maxShared<-max(shared)
	  offsets<-min(as.numeric(grids))
	  for(i in 1:nrow(shared)) {
	    for(j in 1:ncol(shared)) {
		  val<-as.integer(nCol*shared[i,j]/maxShared)
		  xx<-as.numeric(grids[j])
		  yy<-as.numeric(grids[i])
		  rect(xleft=xx-offsets, xright=xx+offsets, ytop=yy+offsets, ybottom=yy-offsets, col=cols[val], border=NA)
		}
	  }

	  colfunc <- colorRampPalette(c("#FFFFFF", "lightblue"))
	  nCol<-100	  
	  cols<-colfunc(nCol)
	  maxD<-max(dronly[,1])
	  offsets<-min(as.numeric(grids))
	  for(i in 1:nrow(dronly)) {
		val<-as.integer(nCol*dronly[i,1]/maxD)
		xx<-as.numeric(grids[i])
		rect(xleft=xx-offsets, xright=xx+offsets, ytop=0-mini_offset, ybottom=xyoffset-mini_offset, col=cols[val], border=NA)
	  }

	  #colfunc <- colorRampPalette(c("#FFFFFF", "papayawhip"))
	  colfunc <- colorRampPalette(c("#FFFFFF", "peachpuff"))
	  nCol<-100	  
	  cols<-colfunc(nCol)
	  maxR<-max(dronly[,2])
	  offsets<-min(as.numeric(grids))
	  for(i in 1:nrow(dronly)) {
		val<-as.integer(nCol*dronly[i,2]/maxR)
		yy<-as.numeric(grids[i])
		rect(xleft=xyoffset-mini_offset, xright=0-mini_offset, ytop=yy+offsets, ybottom=yy-offsets, col=cols[val], border=NA)
	  }
	  
	}
  }
  lines(c(0,1),c(0,0),lty=1,col="#000000")
  lines(c(0,0),c(0,1),lty=1,col="#000000")
  lines(c(0,1),c(1,1),lty=1,col="#000000")
  lines(c(1,1),c(0,1),lty=1,col="#000000")
  axis(1,at=c(0,0.2,0.4,0.6,0.8,1),las=1,cex.axis=0.6)
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),las=1,cex.axis=0.6)
  
  lines(c(0,1),c(0.5,0.5),lty=2,col="#CCCCCC")
  lines(c(0.5,0.5),c(0,1),lty=2,col="#CCCCCC")
  lines(c(xyoffset-mini_offset,0-mini_offset),c(0,0),col="#CCCCCC")
  lines(c(xyoffset-mini_offset,0-mini_offset),c(1,1),col="#CCCCCC")
  lines(c(xyoffset-mini_offset,xyoffset-mini_offset),c(0,1),col="#CCCCCC")
  lines(c(0-mini_offset,0-mini_offset),c(0,1),col="#CCCCCC")
  
  lines(c(0,0),c(xyoffset-mini_offset,0-mini_offset),col="#CCCCCC")
  lines(c(1,1),c(xyoffset-mini_offset,0-mini_offset),col="#CCCCCC")
  lines(c(0,1),c(xyoffset-mini_offset,xyoffset-mini_offset),col="#CCCCCC")
  lines(c(0,1),c(0-mini_offset,0-mini_offset),col="#CCCCCC")
  text(1.025,0-mini_offset-0.015,maxCvg,pos=2,cex=0.5)
  text(1.025,xyoffset-mini_offset+0.005,0,pos=2,cex=0.5)
  text(0-mini_offset-0.02,1.025,maxCvg,pos=1,cex=0.5)
  text(xyoffset-mini_offset+0.01,1.025,0,pos=1,cex=0.5)

  
  text(0.5*xyoffset-0.05,-0.1,info,cex=0.6,pos=1)
  legends<-c(0.1,0.25,0.5,1)*maxCvg
  for(i in 1:length(legends)) if(legends[i]<=maxCvg) {   
    yy<-0-0.2-0.7*i/20
    xx<-0-0.2
    scalee<-maxLen*legends[i]/maxCvg	
    lines(c(xx-scalee,xx),c(yy,yy))
    text(xx,yy,legends[i],pos=4,cex=0.6)
  }
  

  for(i in 1:nrow(onesample)) {
    mind<-as.numeric(onesample[i,"MinD"])
    tind<-as.numeric(onesample[i,"TinD"])
    minr<-as.numeric(onesample[i,"MinR"])
    tinr<-as.numeric(onesample[i,"TinR"])
    PLP<-as.character(onesample[i,"PLP"])
    typee<-as.character(onesample[i,"Type"])
    anno<-""
    if(!(PLP %in% c("P","LP"))) {
      if(onesample[i,"combD"]=="somatic" & onesample[i,"combR"]=="somatic") {
        add_shared(mind,tind,minr,tinr,maxCvg,maxLen,anno,typee)
      }
      if(onesample[i,"combD"]=="wildtype" & onesample[i,"combR"]=="somatic") {
        add_Ronly(mind, tind,minr,tinr,maxCvg,maxLen,leftX=xyoffset-mini_offset,rightX=0-mini_offset,anno,typee)
      }
      if(onesample[i,"combD"]=="somatic" & onesample[i,"combR"]=="wildtype") {
        add_Donly(mind, tind,minr,tinr,maxCvg,maxLen,bottomY=xyoffset-mini_offset,topY=0-mini_offset,anno,typee)
      }
    }
  }
  for(i in 1:nrow(onesample)) {
    mind<-as.numeric(onesample[i,"MinD"])
    tind<-as.numeric(onesample[i,"TinD"])
    minr<-as.numeric(onesample[i,"MinR"])
    tinr<-as.numeric(onesample[i,"TinR"])
    PLP<-as.character(onesample[i,"PLP"])
    typee<-as.character(onesample[i,"Type"])
    anno<-""
    if(PLP %in% c("P","LP")) {
      anno<-paste(c(as.character(onesample[i,"Gene"]),as.character(onesample[i,"AAChange"])),collapse=".")
      if(onesample[i,"combD"]=="somatic" & onesample[i,"combR"]=="somatic") {
        add_shared(mind,tind,minr,tinr,maxCvg,maxLen,anno,typee)
      }
      if(onesample[i,"combD"]=="wildtype" & onesample[i,"combR"]=="somatic") {
        add_Ronly(mind, tind,minr,tinr,maxCvg,maxLen,leftX=xyoffset-mini_offset,rightX=0-mini_offset,anno,typee)
      }
      if(onesample[i,"combD"]=="somatic" & onesample[i,"combR"]=="wildtype") {
        add_Donly(mind, tind,minr,tinr,maxCvg,maxLen,bottomY=xyoffset-mini_offset,topY=0-mini_offset,anno,typee)
      }
    }
  }  
}


data_conversion<-function(ifn.snvinde) {
  dt<-read.table(ifn.snvinde,header=TRUE,sep="\t",colClass="character",check.names=F)
  chr<-matrix(unlist(strsplit(dt[,"Index"],split=".",fixed=TRUE)),ncol=4,byrow=TRUE)[,1]
  pos<-matrix(unlist(strsplit(dt[,"Index"],split=".",fixed=TRUE)),ncol=4,byrow=TRUE)[,2]
  ref<-matrix(unlist(strsplit(dt[,"Index"],split=".",fixed=TRUE)),ncol=4,byrow=TRUE)[,3]
  mut<-matrix(unlist(strsplit(dt[,"Index"],split=".",fixed=TRUE)),ncol=4,byrow=TRUE)[,4]
  SJID<-dt[,"patient"]
  PLP<-dt[,"PLP"]
  Type<-dt[,"mut_type"]

  
  gene<-rep("",nrow(dt))
  clas<-rep("",nrow(dt))
  chag<-rep("",nrow(dt))
  
  MinD1<-rep(0,nrow(dt))
  TinD1<-rep(0,nrow(dt))
  MinR1<-rep(0,nrow(dt))
  TinR1<-rep(0,nrow(dt))
  combD<-rep("wildtype",nrow(dt))
  combR<-rep("wildtype",nrow(dt))

  for(i in 1:nrow(dt)) {
    kk<-unlist(strsplit(as.character(dt[i,"Anno"]),split="\\|\\|"))
	if(length(kk)>=3) {
	  gene[i]<-kk[1]
	  clas[i]<-kk[2]
	  chag[i]<-kk[3]
	}
    if(dt[i,"WGS_D1_manual_call"]=="somatic" & !(dt[i,"WGS_D1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WGS_D1_count"],split=";")))
	  MinD1[i]<-MinD1[i]+zz[1]
	  TinD1[i]<-TinD1[i]+zz[2]
	  combD[i]<-"somatic"
    }
    if(dt[i,"WGS_R1_manual_call"]=="somatic" & !(dt[i,"WGS_R1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WGS_R1_count"],split=";")))
	  MinR1[i]<-MinR1[i]+zz[1]
	  TinR1[i]<-TinR1[i]+zz[2]
	  combR[i]<-"somatic"
    }
    if(dt[i,"WES_D1_manual_call"]=="somatic" & !(dt[i,"WES_D1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WES_D1_count"],split=";")))
	  MinD1[i]<-MinD1[i]+zz[1]
	  TinD1[i]<-TinD1[i]+zz[2]
	  combD[i]<-"somatic"
    }
    if(dt[i,"WES_R1_manual_call"]=="somatic" & !(dt[i,"WES_R1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WES_R1_count"],split=";")))
	  MinR1[i]<-MinR1[i]+zz[1]
	  TinR1[i]<-TinR1[i]+zz[2]
	  combR[i]<-"somatic"
    }

    if(dt[i,"WGS_D1_manual_call"]=="wildtype" & !(dt[i,"WGS_D1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WGS_D1_count"],split=";")))
	  MinD1[i]<-MinD1[i]+0
	  TinD1[i]<-TinD1[i]+zz[2]
    }
    if(dt[i,"WGS_R1_manual_call"]=="wildtype" & !(dt[i,"WGS_R1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WGS_R1_count"],split=";")))
	  MinR1[i]<-MinR1[i]+0
	  TinR1[i]<-TinR1[i]+zz[2]
    }
    if(dt[i,"WES_D1_manual_call"]=="wildtype" & !(dt[i,"WES_D1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WES_D1_count"],split=";")))
	  MinD1[i]<-MinD1[i]+0
	  TinD1[i]<-TinD1[i]+zz[2]
    }
    if(dt[i,"WES_R1_manual_call"]=="wildtype" & !(dt[i,"WES_R1_count"] %in% c("NS"))) {
      zz<-as.integer(unlist(strsplit(dt[i,"WES_R1_count"],split=";")))
	  MinR1[i]<-MinR1[i]+0
	  TinR1[i]<-TinR1[i]+zz[2]
    }  
  }
  ndt<-cbind(SJID,Type=Type,Chr=chr,Pos=pos,Ref=ref,Mut=mut,Gene=gene,Class=clas,AAChange=chag,MinD=MinD1,TinD=TinD1,MinR=MinR1,TinR=TinR1,PLP=PLP,combD,combR)
  #nmss<-setdiff(colnames(ndt),c("quaD1.wes","quaR1.wes","quaD2.wes","quaR2.wes","quaD1.wgs","quaR1.wgs","quaD2.wgs","quaR2.wgs"))
  #nmss<-setdiff(nmss,c("vD1.wes","vR1.wes","vD2.wes","vR2.wes","vD1.wgs","vR1.wgs","vD2.wgs","vR2.wgs","vG1.wes","vG1.wgs","vD1.rna","vR1.rna","vD2.rna","vR2.rna"))
  #nmss<-setdiff(nmss,c("call_D2.wes","call_R2.wes","call_D2.wgs","call_R2.wgs","Chr","Pos"))
  #nmss<-setdiff(nmss,c("cntD2.wes","cntR2.wes","cntD2.wgs","cntR2.wgs","cntG1.wes","cntG1.wgs"))
  #dt<-ndt[,nmss]
  ofn<-gsub(".txt_","_",paste(ifn.snvinde,"_slim.txt",sep=""))
  write.table(ndt,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
  ndt<-read.table(ofn,header=TRUE,sep="\t")
  return(ndt)
}


#setwd("T:/ZhangLab/Projects/COGRelapse/sysRGS01/Oct2/maf2d")
setwd("R:/xmagrp/AMLRelapse/common/COGRelapse/Mar2020/maf2d")

subtypes<-read.table("subtype.txt",header=TRUE,sep="\t",colClass="character")
sjids<-unique(subtypes[,"SJID"])
cloudfiles<-paste("./cloud/",grep("^SJ",dir(path="./cloud",pattern="cloud.txt"),value=TRUE),sep="")

dt<-data_conversion(ifn.snvinde="COGR_coding_SNVIndel.txt")
sv<-read.sv()
burden<-read.table("COGR_SNV_burden.txt",header=TRUE,sep="\t",colClass="character")

pdf("COGRelapse_March2020.pdf",useDingbats=FALSE,width=6,height=6)
par(mar=c(2.5,2.5,0.001,0.001))
for(i in 1:length(sjids)) {
  print(sjids[i])
  if(sjids[i] %in% as.character(dt[,1])) {
    onesample<-dt[as.character(dt[,1])==sjids[i],]
    usi<-as.character(subtypes[subtypes[,"SJID"]==sjids[i],"TID"])
	subtype<-as.character(subtypes[subtypes[,"SJID"]==sjids[i],"March.Major"])
	cloudfn<-grep(sjids[i],cloudfiles,value=TRUE)
    maf2d_one(onesample,sjid=sjids[i],cloudfn=cloudfn,info=paste(c(sjids[i],usi),collapse="/"))
	if(1) {
      if(sjids[i] %in% as.character(sv[,"SJID"])) {
	    shared<-unlist(strsplit(sv[as.character(sv[,"SJID"])==sjids[i],"DR_gene"],split=";"))
	    shared<-c(subtype,shared)
	    Donly<-unlist(strsplit(sv[as.character(sv[,"SJID"])==sjids[i],"D_gene"],split=";"))
	    Ronly<-unlist(strsplit(sv[as.character(sv[,"SJID"])==sjids[i],"R_gene"],split=";"))	
	
	    if(length(shared)>0) {	 text(0.9,0.97,paste(shared, collapse="\n"),cex=0.5,pos=1)	}
	    if(length(Donly)>0) {   text(0.9,-0.05,paste(Donly, collapse="\n"),cex=0.5,pos=1)	}
	    if(length(Ronly)>0) {   text(-0.2,0.97,paste(Ronly, collapse="\n"),cex=0.5,pos=1)	}
	  }
    }
	if(1) {
	  if(sjids[i] %in% as.character(burden[,"SJID"])) {
	    shared<-paste(burden[burden[,"SJID"]==sjids[i],c("nCdsD1R1","nD1R1")],collapse="/")
		Diag<-paste(burden[burden[,"SJID"]==sjids[i],c("nCdsD1","nD1")],collapse="/")
		Relap<-paste(burden[burden[,"SJID"]==sjids[i],c("nCdsR1","nR1")],collapse="/")
		text(0.85,1.02,shared,cex=0.5,pos=4)
		text(0.7,-0.3,Diag,cex=0.5,pos=4)
		text(-0.3,1.02,Relap,cex=0.5,pos=4)
	  }
	}
  }
}
dev.off()

