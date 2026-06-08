

check_it<-function(bamifn,ofn) {
  bams<-readLines(bamifn)
  rlt<-NULL
  for(i in 1:length(bams)) {
    bam <- bams[i]
#    print(bam)
    Sample<-unlist(strsplit(bam, split ="/"))[11]
    Plat<-gsub("TRANSCRIPTOME","RNA",gsub("EXOME","WES",gsub("WHOLE_GENOME","WGS", unlist(strsplit(bam, split ="/"))[12])))  
    SJID<-unlist(strsplit(Sample, split ="_"))[1]
    Smp<-unlist(strsplit(Sample, split ="_"))[2]
    Date<-gsub("-","",unlist(strsplit(as.character(file.info(bam)$mtime),split=" "))[1])
    Year<-unlist(strsplit(unlist(strsplit(as.character(file.info(bam)$mtime),split=" "))[1],split="-"))[1]
    ssize<-as.numeric(file.info(bam)$size)
    if(ssize>10^9) {
      Size<-paste(as.integer(ssize/10^9),"G",sep="")
    } else if(ssize>10^6) {
      Size<-paste(as.integer(ssize/10^6),"M",sep="")
    } else {
      Size<-paste(as.integer(ssize/10^3),"K",sep="")
    }
    codingcvg<-gsub("bam","coverage-post",gsub("\\.bam",".codingexon.cumulative.txt",bam))
    Cvg20X<-"-1"
    Cvg30X<-"-1"
    Cvg40X<-"-1"
    if(file.exists(codingcvg)) {
      if(file.info(codingcvg)$size>0) {
         codingcvgdt<-read.table(codingcvg,header=TRUE,sep="\t",check.names=FALSE)
	 if("19" %in% colnames(codingcvgdt)) {
           Cvg20X<-codingcvgdt["rcumsum","19"]
	 }
	 if("29" %in% colnames(codingcvgdt)) {
           Cvg30X<-codingcvgdt["rcumsum","29"]
	 }
         if("39" %in% colnames(codingcvgdt)) {
           Cvg40X<-codingcvgdt["rcumsum","39"]
         }
      }
    }
    tmp<-c(Sample,bam,Plat,SJID,Smp,Size,Date,Year,Cvg20X,Cvg30X,Cvg40X)
    rlt<-rbind(rlt, tmp)
  }
  colnames(rlt)<-c("Sample","bam","Plat","SJID","Smp","Size","Date","Year","Cvg20X","Cvg30X","Cvg40X")
  write.table(rlt, ofn, row.names=FALSE,col.names=TRUE,sep="\t",quote=F)
}


args<-commandArgs(trailingOnly=TRUE)
if(length(args)==2) {
  bamifn<-args[1]
  ofn<-args[2]
  check_it(bamifn,ofn)
} else {
  print("Rscript ~ <bamlist.ifn> <info.ofn>")
}

