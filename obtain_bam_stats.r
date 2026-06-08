obtain_machine_flowcellname<-function(bam,prefix) {
  tmpreadnamefile<-paste(c(prefix,"_tmp_readnames.txt"),collapse="")
  cmd<-paste(c("samtools view ", bam, " | head -n10000 | cut -f1 > ", tmpreadnamefile),collapse=" ")
  system(cmd)
  lns<-readLines(tmpreadnamefile)
  dt<-matrix(unlist(strsplit(lns,split=":")),nrow=length(lns),byrow=TRUE)
  machine<-"invalid_machine"
  flowcell<-"invalid_readname"
  if(ncol(dt)==5) {
    flowcell<-paste(unique(dt[,1]),collapse=";")
  }
  if(ncol(dt)==7) {
    machine<-paste(sort(unique(dt[,1])),collapse=";")
    flowcell<-paste(sort(unique(dt[,3])),collapse=";")
  }
  ret<-c(flowcell,machine)
  system(paste(c("rm ",tmpreadnamefile),collapse=" "))
  return(ret)
}

obtain_library_name<-function(bam,prefix) {
  tmpreadnamefile<-paste(c(prefix,"_tmp_libnames.txt"),collapse="")
  cmd<-paste(c("samtools view -H ", bam, " | grep -e SM: -e SN:1 -e SN:chr1 > ",tmpreadnamefile),collapse=" ")
  system(cmd)
  dt<-readLines(tmpreadnamefile)
  lib<-NULL
  for(i in 1:length(dt)) {
    lib<-union(lib,gsub("SM:","",grep("SM:",unlist(strsplit(dt[i],split="\t")),value=TRUE)))
  }
  if(is.null(lib)) {
    ret<-"NoLibraryName"
  } else {
    ret<-paste(lib,collapse=";")
  }
  system(paste(c("rm ",tmpreadnamefile),collapse=" "))
  return(ret)
}

obtain_build<-function(bam,prefix) {
  tmpreadnamefile<-paste(c(prefix,"_tmp_libnames.txt"),collapse="")
  cmd<-paste(c("samtools view -H ", bam, " | grep -e SN:1 -e SN:chr1 | grep LN: > ",tmpreadnamefile),collapse=" ")
  system(cmd)
  dt<-readLines(tmpreadnamefile)
  chr1_len<-NULL
  ret<-"NoBuild"
  for(i in 1:length(dt)) {
    chrr<-grep("SN:",unlist(strsplit(dt[i],split="\t")),value=TRUE)
    if(chrr %in% c("SN:chr1","SN:1")) {
      chr1_len<-union(chr1_len,gsub("LN:","",grep("LN:",unlist(strsplit(dt[i],split="\t")),value=TRUE)))
    }
  }
  if(!is.null(chr1_len)) {
    if(chr1_len=="248956422") ret<-"hg38"
    if(chr1_len=="249250621") ret<-"hg19"
    if(chr1_len=="248387328") ret<-"hgT2T"
    if(chr1_len=="247249719") ret<-"hg18"
  }
  system(paste(c("rm ",tmpreadnamefile),collapse=" "))
  return(ret)
}



check_it<-function(bamlistifn,ofn) {
  bams<-readLines(bamlistifn)
  rlt<-NULL
  for(i in 1:length(bams)) {
    bam <- bams[i]
    print(paste(c(i, length(bams)),collapse= " / "))
    Flowcell_Machine<-obtain_machine_flowcellname(bam,prefix=bamlistifn)
    LibName<-obtain_library_name(bam,prefix=bamlistifn)
    Build<-obtain_build(bam,prefix=bamlistifn)
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
    tmp<-c(Sample,bam,Plat,SJID,Smp,Size,Date,Year,Cvg20X,Cvg30X,Cvg40X,Flowcell_Machine,LibName,Build)
    rlt<-rbind(rlt, tmp)
  }
  colnames(rlt)<-c("Sample","bam","Plat","SJID","Smp","Size","Date","Year","Cvg20X","Cvg30X","Cvg40X","Flowcell","Machine","LibName","Build")
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

